using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using MathNet.Numerics.Integration;

namespace Master.Components
{
    public class KMatrix : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the kMatrix class.
        /// </summary>
        public KMatrix()
          : base("KMatrix", "Nickname",
              "Description",
              "Løve", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("bars","b","bars of the geometry",GH_ParamAccess.list);
            
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("kMatrix","ke","ke created from MathCad",GH_ParamAccess.list);
            pManager.AddMatrixParameter("Nmatrix", "N", "N created from MathCad", GH_ParamAccess.item);
            pManager.AddMatrixParameter("dNMatrix", "dN", "kdN created from MathCad", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input

            List<BarClass> bars = new List<BarClass>();

            DA.SetDataList(0, bars);

            //code

            CreateGenerelizedShapeFunc(bars, out Matrix<double> N, out Matrix<double> dN);

            //Matrix<double> ke = DenseMatrix.OfArray(new double[12, 12]);

            List<Matrix<double>> L_ke = new List<Matrix<double>>();

            foreach (BarClass b in bars)
            {
                Line currentLine = b.axis;
                double L = currentLine.Length * 1000;

                double X = b.section.CSA * b.material.youngsModolus / L;
                double Y1 = 12.00 * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 3));
                double Y2 = 6.00 * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 2));  //slender profile, approx.
                double Y3 = 4.00 * b.material.youngsModolus * b.section.Iz / L;
                double Y4 = 2.00 * b.material.youngsModolus * b.section.Iz / L;

                double Z1 = 12.00 * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 3));
                double Z2 = 6.00 * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 2));  //slender profile, approx.
                double Z3 = 4.00 * b.material.youngsModolus * b.section.Iy / L;
                double Z4 = 2.00 * b.material.youngsModolus * b.section.Iy / L;

                double C = (b.material.G * b.section.J) / L;

                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                {
                {X,  0,   0,   0,   0,   0,  -X,   0,   0,   0,   0,   0},
                {0,  Y1,  0,   0,   0,   Y2,  0,  -Y1,  0,   0,   0,   Y2},
                {0,  0,   Z1,  0,  -Z2,  0,   0,   0,  -Z1,  0,  -Z2,  0},
                {0,  0,   0,   C,   0,   0,   0,   0,   0,  -C,   0,   0},
                {0,  0 , -Z2,  0,   Z3,  0,   0,   0,   Z2,  0,   Z4,  0},
                {0,  Y2,  0,   0,   0,   Y3,  0,  -Y2,  0,   0,   0,   Y4},
                {-X, 0,   0,   0,   0,   0,   X,   0,   0,   0,   0,   0},
                {0, -Y1,  0,   0,   0,  -Y2,  0,   Y1,  0,   0,   0,  -Y2},
                {0,  0,  -Z1,  0,   Z2,  0,   0,   0,   Z1,  0,   Z2,  0},
                {0,  0,   0,  -C,   0,   0,   0,   0,   0,   C,   0,   0},
                {0,  0,  -Z2,  0,   Z4,  0,   0,   0,   Z2,  0,   Z3,  0},
                {0,  Y2,  0,   0,   0,   Y4,  0,  -Y2,  0,   0,   0,   Y3},


                });

                L_ke.Add(ke);
                
            }
            

            DA.SetDataList(0, L_ke);
            DA.SetData(1, N);
            DA.SetData(2, dN);
        }

        private static void CreateGenerelizedShapeFunc(List<BarClass> _bars, out Matrix<double> N, out Matrix<double> dN)
        {

            Matrix<double> _N = DenseMatrix.OfArray(new double[6, 12]);
            Matrix<double> _dN = DenseMatrix.OfArray(new double[6, 12]);
            foreach (BarClass b in _bars)
            {
                Line currentLine = b.axis;
                double L = currentLine.Length;
                var x = L / 2;

                double N1 = 1 - x / L;                                                                      //axial translation in node 1
                double N2 = x / L;                                                                          //axial translation in node 2
                double N3 = 1 - 3 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2 * Math.Pow(x, 3) / Math.Pow(L, 3);  //translation node 1
                double N4 = -x * (1 - 2 * x / L + Math.Pow(x, 2) / Math.Pow(L, 2));                                 //rotation node 1
                double N5 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * Math.Pow(x, 3) / Math.Pow(L, 3);      //translation node 2
                double N6 = x * (x / L - Math.Pow(x, 2) / Math.Pow(L, 2));                                 //rotation node 2


                double dN1 = -1 / L;
                double dN2 = 1 / L;
                double dN3 = -6 * x / Math.Pow(L, 2) + 6 * Math.Pow(x, 2) / Math.Pow(L, 3);
                double dN4 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 4 * x / L + 1;
                double dN5 = -dN3;
                double dN6 = 3 * Math.Pow(x, 2) / Math.Pow(L, 2) - 2 * x / L;

                double ddN3 = -6 / Math.Pow(L, 2) + 12 * x / Math.Pow(L, 3);
                double ddN4 = 6 * x / Math.Pow(L, 2) - 4 / L;
                double ddN5 = 6 / Math.Pow(L, 2) - 12 * x / Math.Pow(L, 3);
                double ddN6 = 6 * x / Math.Pow(L, 2) - 2 / L;


                // første rad blir ux = N1*ux1 + N2*ux2
                //andre rad blir uy = N3*uy1 + N5*uy2 etc.

                _N = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                {N1,   0,     0,    0,    0,     0,    N2,     0,     0,     0,     0,    0},
                {0,   N3,     0,    0,    0,    N4,     0,    N5,     0,     0,     0,   N6},
                {0,    0,    N3,    0,  -N4,     0,     0,     0,    N5,     0,   -N6,    0},
                {0,    0,     0,   N1,    0,     0,     0,     0,     0,    N2,     0,    0},


                {0,    dN3,   0,    0,    0,   dN4,     0,    dN5,    0,     0,     0,   dN6},
                {0,     0,   dN3,   0,  -dN4,    0,     0,     0,    dN5,    0,   -dN6,   0},

                });




                _dN = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                {dN1,   0,    0,    0,    0,     0,    dN2,    0,     0,     0,     0,    0},
                {0,    dN3,   0,    0,    0,   dN4,     0,    dN5,    0,     0,     0,   dN6},
                {0,     0,   dN3,   0,  -dN4,    0,     0,     0,    dN5,    0,   -dN6,   0},
                {0,     0,    0,   dN1,   0,     0,     0,     0,     0,    dN2,    0,    0},
                {0,    ddN3,  0,    0,    0,   ddN4,    0,    ddN5,   0,     0,     0,   ddN6},
                {0,     0,   ddN3,   0,  -ddN4,    0,   0,     0,    ddN5,   0,   -ddN6,   0},

                });



            }

            N = _N;
            dN = _dN;
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
            get { return new Guid("9B084A55-ED63-4E33-A5FC-D0310EBF2545"); }
        }
    }
}