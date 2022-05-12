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
            pManager.AddMatrixParameter("ke", "ke", "Element stiffness matrix for NURBS element", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            List<Point3d> ipts = new List<Point3d>();
            int ngp = new int();
            List<double> X = new List<double>();
            int p = new int();



            DA.GetDataList(0, ipts);
            DA.GetData(1, ref ngp);
            DA.GetDataList(2, X);
            DA.GetData(3, ref p);


            List<double> Xi = new List<double>();
            for (int i=0; i<ngp; i++)
            {
                Xi.Add(i / (ngp - 1));
            }

            foreach(double xi in Xi)
            {
                //info
                var info = new List<string>() { "code for NURBS" };

                //control points

                var cpts = ipts;


                //knot vector
                List<double> kV = X;// new List<double>(){0,0,0,0.5,1,1,1};

                string kVstring = "KnotVec = ";
                foreach (var v in kV)
                    kVstring = kVstring + " " + v;
                info.Add(kVstring);

                //degree
                int pu = 1;
                int pv = p;
                int n = cpts.Count;

                //double xi = 0.49;

                //basic shape function degree 0
                List<double> N0 = getN0(kV, xi);
                string N0string = "N0 = ";
                foreach (var n0 in N0)
                    N0string = N0string + " " + n0;
                info.Add(N0string);

                //basic shape function degree 1 and higher
                int d = 1;
                List<double> N = getN(kV, xi, p);
                List<double> dN = getdN(kV, xi, p, d);
                string Nstring = "N = ";
                foreach (var nshp in N)
                    Nstring = Nstring + " " + nshp;
                info.Add(Nstring);

                //Beam properties
                double L = 1.3;
                double J = L / 2;
                double A = 0.2;
                double I = 0.014;
                double E = 21000;


                //Selective basis functions (for vertical disp and rot)

                int db = 2;   //derived d-times
                double B0 = getdN(kV, xi, p, d)[0] + getdN(kV, xi, p, db)[1];
                double B1 = getdN(kV, xi, p, d)[2] + getdN(kV, xi, p, db)[3];
                double B2 = L / 3 * getdN(kV, xi, p, db)[1];
                double B3 = -L / 3 * getdN(kV, xi, p, db)[2];


                //Making K-matrix


                List<double> GPW = new List<double>() { 2, 1, }; 
                int dim = pu + 1 + pv + 1;
                Matrix<double> k_e = DenseMatrix.OfArray(new double[dim, dim]);

                for (int a = 0; a < N.Count; a++)
                {
                    for (int b = 0; b < N.Count; b++)
                    {
                        k_e.Column(a)[b] += (E * A / J * J) * dN[a] * dN[b] * J * GPW[ngp - 1];
                    }
                }
                



    //calculate
                int ispan = findSpan(kV, xi);
                info.Add("ispan for xi " + xi + " is =" + ispan);

                double posX = 0;
                double posY = 0;
                double posZ = 0;

                for (int c = 0; c < cpts.Count; c++)
                {
                    posX = posX + cpts[c].X * N[c];
                    posY = posY + cpts[c].Y * N[c];
                    posZ = posZ + cpts[c].Z * N[c];
                }

                Point3d pos = new Point3d(posX, posY, posZ);
                info.Add("Xcord = " + posX);
                info.Add("Ycord = " + posX);
                info.Add("Zcord = " + posX);


                infos = info;
                opts = cpts;
                opt = pos;
            }
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

            if (d == 1)
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
            else if (p > 1)
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