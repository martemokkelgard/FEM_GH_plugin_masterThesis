using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BeamClass2D
    {
        public string Name;
        public int Id;
        public SectionClass2D section;
        public MaterialClass2D material;
        public Line axis;
        public NodeClass2D startNode;
        public NodeClass2D endNode;



        public BeamClass2D()
        { }

        public BeamClass2D(string _Name, Line _line, SectionClass2D _section, MaterialClass2D _material)
        {
            Name = _Name;
            axis = _line;
            section = _section;
            material = _material;
        }
    }
}