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
    internal class BeamClassNonLin
    {
        public string Name;
        public int Id;
        public SectionClass section;
        public MaterialClass material;
        public Curve axis;
        public List<NodeClass> nodes;
       


        public BeamClassNonLin()
        { }

        public BeamClassNonLin(string _Name, Curve _line, SectionClass _section, MaterialClass _material)
        {
            Name = _Name;
            axis = _line;
            section = _section;
            material = _material;
        }
    }
}
