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
using MathNet.Symbolics;


namespace Master.Components.Nonlinear
{
    public class K_matrix : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the K_matrix class.
        /// </summary>
        public K_matrix()
          : base("K-matrix", "Nickname",
              "Description",
              "Panda", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Degree","deg","Degree of nonlinearity",GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("MatrixClass", "mC", "Matrix class with degree and element matrix", GH_ParamAccess.item);
           
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            /*
            
            int deg = 1;
            int dofs = 6 * (deg + 1);


            DA.GetData(0, ref deg);

            Vector<double> Nq = DenseVector.OfArray(new double[dofs]);
            var x = SymbolicExpression.Variable("x");
            Convert.ToDouble(x);


            for (int i = 0; i < (deg+1) *2; i++)
            {
                Nq[i] = Math.Pow(x, i);
            }

            

            Matrix<double> ke = DenseMatrix.OfArray(new double[dofs, dofs]);    //empty ke stiffness matrix

            for (int i = 0; i < (deg +1); i++)
            {
                for (int j = 0; j < (deg +1); j++)
                {
                    ke[i * 6, j * 6] = k1[i, j];
                    ke[i * 6 + 3, j * 6 + 3] = k2[i, j];

                    ke[i * 6 + 1, j * 6 + 1] = ky[i * 2, j * 2];
                    ke[i * 6 + 1, j * 6 + 5] = ky[i * 2, j * 2 + 1];

                    ke[i * 6 + 2, j * 6 + 2] = kz[i * 2, j * 2];
                    ke[i * 6 + 2, j * 6 + 4] = kz[i * 2, j * 2 + 1];

                    ke[i * 6 + 4, j * 6 + 2] = kz[i * 2 + 1, j * 2];
                    ke[i * 6 + 4, j * 6 + 4] = kz[i * 2 + 1, j * 2 + 1];

                    ke[i * 6 + 5, j * 6 + 1] = ky[i * 2 + 1, j * 2];
                    ke[i * 6 + 5, j * 6 + 5] = ky[i * 2 + 1, j * 2 + 1];
                }
            }
            */
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
            get { return new Guid("7A182C6B-AFC8-4803-A860-635DDAEAD98D"); }
        }
    }
}