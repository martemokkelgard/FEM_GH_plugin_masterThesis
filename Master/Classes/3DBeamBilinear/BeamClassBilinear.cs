using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BeamClassBilinear
    {
        public string Name;
        public int Id;
        public SectionClass section;
        public MaterialClass material;
        public Curve axis;
        public NodeClass startNode;
        public NodeClass endNode;
        public NodeClass midNode;



        public BeamClassBilinear()
        { }

        public BeamClassBilinear(string _Name, Curve _line, SectionClass _section, MaterialClass _material)
        {
            Name = _Name;
            axis = _line;
            section = _section;
            material = _material;
        }
    }
}
