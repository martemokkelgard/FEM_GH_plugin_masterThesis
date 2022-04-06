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

namespace Master
{
    internal class MatrixClass
    {
        public Matrix<double> elementMatrix;
        public double deg;


        public MatrixClass()
        {

        }

        public MatrixClass(double _deg, Matrix<double> _elementMatrix)
        {
            elementMatrix = _elementMatrix;
            deg = _deg;
        }
    }
}
