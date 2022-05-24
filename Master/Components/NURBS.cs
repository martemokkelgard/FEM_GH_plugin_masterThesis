using System;
using System.Collections.Generic;
using System.Windows;
using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using MathNet.Numerics.Integration;
using System.Windows.Media;
namespace Master.Components
{
    public class NURBS : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the NURBS class.
        /// </summary>
        public NURBS()
          : base("NURBS", "Nickname",
              "Description",
              "Panda", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("ControlPoints", "CP", "Control Points", GH_ParamAccess.list);
            pManager.AddNumberParameter("ngp", "ngp", "Number of Gauss-points", GH_ParamAccess.item);
            pManager.AddNumberParameter("Knot", "Knot", "Knot Vector", GH_ParamAccess.list);
            pManager.AddNumberParameter("Deg", "Deg", "Degree", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("ke", "ke", "Element stiffness matrix for NURBS element", GH_ParamAccess.item);
            pManager.AddPointParameter("CurvePoint", "CurvePoint", "Point on Curve", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            List<Point3d> ipts = new List<Point3d>();
            double ngp = new double();
            List<double> X = new List<double>();
            double pp = new double();



            DA.GetDataList(0, ipts);
            DA.GetData(1, ref ngp);
            DA.GetDataList(2, X);
            DA.GetData(3, ref pp);

            ngp = Convert.ToInt32(ngp);
            int p = Convert.ToInt32(pp);

            List<Point3d> position = new List<Point3d>();
            //List<double> GPW = new List<double>();

            var GPW = getGPW(Convert.ToInt32(ngp));


            /*
            if (ngp == 2)
            {
                Xi.Add(-Math.Sqrt(1.0 / 3.0));
                Xi.Add(Math.Sqrt(1.0 / 3.0));
                GPW.Add(1.0);
                GPW.Add(1.0);
            }
            */

            //control points

            var cpts = ipts;

            // degree
            int pu = 1;
            int pv = p;
            int cpc = cpts.Count;

            int dim = pu + 1 + pv + 1;
            Matrix<double> k_e = DenseMatrix.OfArray(new double[dim, dim]);
            Matrix<double> k_ee = DenseMatrix.OfArray(new double[4, 4]);

            for (int n = 0; n < ngp; n++)
            {
                //info
                var info = new List<string>() { "code for NURBS" };


                //knot vector
                List<double> kVv = X;// new List<double>(){0,0,0,0.5,1,1,1};
                List<double> kVu = new List<double>() { -1, -1, 1, 1 };

                string kVstring = "KnotVec = ";
                foreach (var v in kVv)
                    kVstring = kVstring + " " + v;
                info.Add(kVstring);


                //double xi = 0.49;

                //basic shape function degree 0
                List<double> N0v = getN0(kVv, GPW[0, n]);
                List<double> N0u = getN0(kVu, GPW[0, n]);
                string N0string = "N0 = ";
                foreach (var n0 in N0v)
                    N0string = N0string + " " + n0;
                info.Add(N0string);

                //basic shape function degree 1 and higher
                int d = 1;
                List<double> Nu = getN(kVu, GPW[0, n], pu);
                List<double> dNu = getdN(kVu, GPW[0, n], pu, d);
                List<double> Nv = getN(kVv, GPW[0, n], pv);
                List<double> dNv = getdN(kVv, GPW[0, n], pv, d);
                string Nstring = "N = ";
                foreach (var nshp in Nv)
                    Nstring = Nstring + " " + nshp;
                info.Add(Nstring);

                //Beam properties
                double L = 1.3;
                double J = L / 2;
                double A = 0.2;
                double I = 0.014;
                double E = 1525;

                /*
                //Selective basis functions (for vertical disp and rot)

                int db = 2;   //derived d-times
                List<double> B = new List<double>();
                double B0 = getdN(kVv, GPW[0, n], pv, db)[0] + getdN(kVv, GPW[0, n], pv, db)[1];
                double B1 = getdN(kVv, GPW[0, n], pv, db)[2] + getdN(kVv, GPW[0, n], pv, db)[3];
                double B2 = L / 3.0 * getdN(kVv, GPW[0, n], pv, db)[1];
                double B3 = -L / 3.0 * getdN(kVv, GPW[0, n], pv, db)[2];
                B.Add(B0); B.Add(B1); B.Add(B2); B.Add(B3);


                //Making K-matrix

                


                for (int a = 0; a < Nu.Count; a++)
                {
                    for (int b = 0; b < Nu.Count; b++)
                    {
                        k_e[3 * a, 3 * b] += (E * A / (J * J)) * dNu[a] * dNu[b] * J * GPW[1, n];
                    }
                }


                for (int g = 0; g < Nv.Count; g++)
                {
                    for (int e = 0; e < Nv.Count; e++)

                    {
                        k_ee[g, e] += (E * I) / (Math.Pow(J, 4)) * B[g] * B[e] * J * GPW[1, n];

                    }
                }


                */

                //calculate
                int ispan = findSpan(kVv, GPW[0, n]);
                info.Add("ispan for xi " + GPW[0, n] + " is =" + ispan);

                double posX = 0;
                double posY = 0;
                double posZ = 0;

                for (int c = 0; c < cpts.Count; c++)
                {
                    posX = posX + cpts[c].X * Nv[c];
                    posY = posY + cpts[c].Y * Nv[c];
                    posZ = posZ + cpts[c].Z * Nv[c];
                }

                Point3d pos = new Point3d(posX, posY, posZ);
                info.Add("Xcord = " + posX);
                info.Add("Ycord = " + posX);
                info.Add("Zcord = " + posX);

                position.Add(pos);

            }

            //Curve crv = Curve.CreateInterpolatedCurve(position, p);

            DA.SetData(0, k_ee);
            DA.SetDataList(1, position);
        }







        private int findSpan(List<double> kV, double xi)
        {
            int ispan = 0;

            for (int i = 0; i < kV.Count - 1; i++)
            {
                double span0 = kV[i];
                double span1 = kV[i + 1];

                if (xi >= span0 && xi < span1)
                { ispan = i; }
                else if (xi > span0 && xi <= span1)
                { ispan = i; }
            }
            return ispan;
        }

        private double[,] getGPW(int ngp)
        {
            double[,] gpw = new double[,] { };
            switch (ngp)
            {
                case 1:
                    double p11 = 0;
                    double w11 = 2;
                    gpw = new double[,] { { p11 }, { w11 } };
                    break;
                case 2:
                    double p21 = -Math.Sqrt(0.33);
                    double w21 = 1;
                    double p22 = Math.Sqrt(0.33);
                    double w22 = 1;
                    gpw = new double[,] { { p21, p22 }, { w21, w22 } };
                    break;
                case 3:
                    double p31 = -Math.Sqrt(0.6);
                    double w31 = 0.55556;
                    double p32 = 0;
                    double w32 = 0.88889;
                    double p33 = Math.Sqrt(0.6);
                    double w33 = 0.55556;
                    gpw = new double[,] { { p31, p32, p33 }, { w31, w32, w33 } };
                    break;
                case 4:
                    double p41 = -Math.Sqrt(0.42857 + 0.28571 * Math.Sqrt(1.2));
                    double w41 = (18 - Math.Sqrt(30)) / 36.0;
                    double p42 = -Math.Sqrt(0.42857 - 0.28571 * Math.Sqrt(1.2));
                    double w42 = (18 + Math.Sqrt(30)) / 36.0;
                    double p43 = Math.Sqrt(0.42857 - 0.28571 * Math.Sqrt(1.2));
                    double w43 = (18 + Math.Sqrt(30)) / 36.0;
                    double p44 = Math.Sqrt(0.42857 + 0.28571 * Math.Sqrt(1.2));
                    double w44 = (18 - Math.Sqrt(30)) / 36.0;
                    gpw = new double[,] { { p41, p42, p43, p44 }, { w41, w42, w43, w44 } };
                    break;
                case 5:
                    double p51 = -0.33333 * Math.Sqrt(5 + 2 * Math.Sqrt(1.42857));
                    double w51 = (322 - 13 * Math.Sqrt(70)) / 900;
                    double p52 = -0.33333 * Math.Sqrt(5 - 2 * Math.Sqrt(1.42857));
                    double w52 = (322 + 13 * Math.Sqrt(70)) / 900;
                    double p53 = 0;
                    double w53 = (128 / 225.0);
                    double p54 = 0.33333 * Math.Sqrt(5 - 2 * Math.Sqrt(1.42857));
                    double w54 = (322 + 13 * Math.Sqrt(70)) / 900;
                    double p55 = 0.33333 * Math.Sqrt(5 + 2 * Math.Sqrt(1.42857));
                    double w55 = (322 - 13 * Math.Sqrt(70)) / 900;
                    gpw = new double[,] { { p51, p52, p53, p54, p55 }, { w51, w52, w53, w54, w55 } };
                    break;



            }
            return gpw;



        }
        private List<double> getN(List<double> kV, double xi, int p)
        {
            List<double> N = new List<double>();

            List<double> N0 = getN0(kV, xi); //get basic shape function

            if (p == 1)
            {
                for (int i = 0; i < kV.Count - p - 1; i++)
                {
                    double a1 = xi - kV[i];
                    double a2 = kV[i + p] - kV[i];
                    double a = 0;
                    if (a2 != 0)
                        a = a1 / a2;

                    double b1 = kV[i + p + 1] - xi;
                    double b2 = kV[i + p + 1] - kV[i + 1];
                    double b = 0;
                    if (b2 != 0)
                        b = b1 / b2;

                    double n = a * N0[i] + b * N0[i + 1];
                    N.Add(n);

                }
            }
            else if (p > 1)
            {
                for (int i = 0; i < kV.Count - p - 1; i++)
                {
                    double a1 = xi - kV[i];
                    double a2 = kV[i + p] - kV[i];
                    double a = 0;
                    if (a2 != 0)
                        a = a1 / a2;

                    double b1 = kV[i + p + 1] - xi;
                    double b2 = kV[i + p + 1] - kV[i + 1];
                    double b = 0;
                    if (b2 != 0)
                        b = b1 / b2;

                    double n = a * getN(kV, xi, p - 1)[i] + b * getN(kV, xi, p - 1)[i + 1];
                    N.Add(n);
                }
            }


            return N;
        }


        private List<double> getdN(List<double> kV, double xi, int p, int d)
        {
            List<double> dN = new List<double>();

            List<double> N0 = getN0(kV, xi); //get basic shape function


            if (d == 1 && p == 1)
            {
                for (int i = 0; i < kV.Count - p - 1; i++)
                {
                    double a1 = p;
                    double a2 = kV[i + p] - kV[i];
                    double a = 0;
                    if (a2 != 0)
                        a = a1 / a2;

                    double b1 = p;
                    double b2 = kV[i + p + 1] - kV[i + 1];
                    double b = 0;
                    if (b2 != 0)
                        b = b1 / b2;

                    double n = a * N0[i] + b * N0[i + 1];
                    dN.Add(n);

                }
            }

            else if (d == 1 && p > 1)
            {
                for (int i = 0; i < kV.Count - p - 1; i++)
                {
                    double a1 = p;
                    double a2 = kV[i + p] - kV[i];
                    double a = 0;
                    if (a2 != 0)
                        a = a1 / a2;

                    double b1 = p;
                    double b2 = kV[i + p + 1] - kV[i + 1];
                    double b = 0;
                    if (b2 != 0)
                        b = b1 / b2;

                    double n = a * getN(kV, xi, p - 1)[i] + b * getN(kV, xi, p - 1)[i + 1];
                    dN.Add(n);

                }
            }
            else if (d > 1 && p > 1)
            {
                for (int i = 0; i < kV.Count - p - 1; i++)
                {
                    double a1 = p;
                    double a2 = kV[i + p] - kV[i];
                    double a = 0;
                    if (a2 != 0)
                        a = a1 / a2;

                    double b1 = p;
                    double b2 = kV[i + p + 1] - kV[i + 1];
                    double b = 0;
                    if (b2 != 0)
                        b = b1 / b2;

                    double n = a * getdN(kV, xi, p - 1, d - 1)[i] + b * getdN(kV, xi, p - 1, d - 1)[i + 1];
                    dN.Add(n);
                }
            }


            return dN;
        }


        private List<double> getN0(List<double> kV, double xi)
        {
            List<double> N0 = new List<double>();

            for (int i = 0; i < kV.Count - 1; i++)
            {
                double xi0 = kV[i];     //xi for i
                double xi1 = kV[i + 1]; //xi for i+1
                if (xi >= xi0 && xi < xi1)
                { N0.Add(1); }
                else if (xi == 1 && xi1 == 1 && xi0 != 1)
                { N0.Add(1); }
                else
                { N0.Add(0); }
            }
            return N0;
        }





        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2E4586E7-E63A-43EA-A741-ABAAB45F420F"); }
        }
    }
}