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
            pManager.AddNumberParameter("Knot axial", "Knot", "Knot Vector", GH_ParamAccess.list);
            pManager.AddNumberParameter("Knot vert and rot", "Knot", "Knot Vector", GH_ParamAccess.list);
            pManager.AddMatrixParameter("Deg", "Deg", "Degree", GH_ParamAccess.item);
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

            List<Point3d> ipts = new List<Point3d>();   //control points
            double ngpp = new double();                        //Number of gauss points
            List<double> X_ax = new List<double>();        //Knot vector
            List<double> X_ver = new List<double>();
            double pp = new double();                          //degree



            DA.GetDataList(0, ipts);
            DA.GetData(1, ref ngpp);
            DA.GetDataList(2, X_ax);
            DA.GetDataList(3, X_ver);
            DA.GetData(4, ref pp);

            int ngp = Convert.ToInt32(ngpp);
            int p = Convert.ToInt32(pp);

            List<double> Xi = new List<double>();
            Xi.Add(-Math.Sqrt(1.0 / 3.0));
            Xi.Add(Math.Sqrt(1.0 / 3.0));


            int pu = 1;
            int pv = 3;
            int dim = pu + 1 + pv + 1;
            Matrix<double> k_e = DenseMatrix.OfArray(new double[dim, dim]);

            foreach (double xi in Xi)
            {
                //info
                var info = new List<string>() { "code for NURBS" };

                //control points

                var cpts = ipts;


                //knot vector
                List<double> kV_ax = X_ax;// new List<double>(){-1,-1,1,1};
                List<double> kV_vert = X_ver;// new List<double>(){-1,-1,-1,-1,1,1,1,1};


                //degree
                int n = cpts.Count;

                //double xi = 0.49;

                //basic shape function degree 0
                Vector<double> N0_ax = getN0(kV_ax, xi);


                //basic shape function degree 1 and higher
                int d = 1;
                Vector<double> N_ax = getN(kV_ax, xi, 1);
                Vector<double> dN_ax = getdN(kV_ax, xi, 1, d);

                Vector<double> dN_vert = getdN(kV_ax, xi, 1, d);
                string Nstring = "N = ";


                //Beam properties
                double L = 1.3;
                double J = L / 2.0;
                double A = 0.2;
                double I = 0.014;
                double E = 1525.0;


                //Selective basis functions (for vertical disp and rot)

                int db = 2;   //derived d-times
                


                //Making K-matrix


                List<double> GPW = new List<double>() { 2, 1, }; 
                

                for (int a = 0; a < ngp; a++)
                {
                    for (int b = 0; b < 2; b++)
                    {
                        for(int c = 0; b < 2; b++)
                        {
                            k_e[3 * b, 3 * c] += (E * A / J * J) * dN_ax[c] * dN_ax[b] * J;
                        }
                            
                    }
                }

                Vector<double> B_vec = SparseVector.OfEnumerable(new double[4]);

                for (int a = 0; a < ngp; a++)
                {
                    B_vec[0] = getdN(kV_vert, xi, pv, d)[0] + getdN(kV_vert, xi, pv, db)[1];
                    B_vec[1] = getdN(kV_vert, xi, pv, d)[2] + getdN(kV_vert, xi, pv, db)[3];
                    B_vec[2] = L / 3 * getdN(kV_vert, xi, p, db)[1];
                    B_vec[3] = -L / 3 * getdN(kV_vert, xi, p, db)[2];

                    for (int j = 1; j < 4; j++)
                    {
                        for (int k = 1; k < 4; k++)
                        {
                            k_e[j+2, k+2] += (E * I * B_vec[j] * B_vec[k] * J);
                        }

                    }
                }



                //calculate
                /*
                            int ispan = findSpan(kV_ax, xi);
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


                            //infos = info;
                            //opts = cpts;
                            //opt = pos;
                */

            }

            DA.SetData(0, k_e);



        }

        private int findSpan(List<double> kV, double xi)
        {
            int ispan = 0;

            for (int i = 0; i < kV.Count -1; i++)
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


        private Vector<double> getN(List<double> kV, double xi, int p)
        {
            Vector<double> N = SparseVector.OfEnumerable(new double[kV.Count-p-1]);

            Vector<double> N0 = getN0(kV, xi); //get basic shape function

            if (p == 1)
            {
                for (int i = 0; i < kV.Count - p-1; i++)
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
                    N[i] = (n);

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
                    N[i] = (n);
                }
            }


            return N;
        }


        private Vector<double> getdN(List<double> kV, double xi, int p, int d)
        {
            Vector<double> dN = SparseVector.OfEnumerable(new double[kV.Count-p-1]);

            Vector<double> N0 = getN0(kV, xi); //get basic shape function

            if (p == 1)
            {
                for (int i = 0; i < kV.Count - p-1; i++)
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
                    dN[i] = (n);

                }
            }
            else if (p > 1)
            {
                for (int i = 0; i < kV.Count -p - 1; i++)
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
                    dN[i] = (n);
                }
            }


            return dN;
        }


        private Vector<double> getN0(List<double> kV, double xi)
        {
            Vector<double> N0 = SparseVector.OfEnumerable(new double[kV.Count-1]);

            for (int i = 0; i < kV.Count - 1; i++)
            {
                double xi0 = kV[i];     //xi for i
                double xi1 = kV[i + 1]; //xi for i+1
                if (xi >= xi0 && xi < xi1)
                { N0[i] = (1); }
                else if (xi == 1 && xi1 == 1 && xi0 != 1)
                { N0[i] = (1); }
                else
                { N0[i] = (0); }
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