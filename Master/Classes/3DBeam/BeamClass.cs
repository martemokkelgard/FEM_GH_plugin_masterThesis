using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BeamClass
    {
        public string Name;
        public int Id;
        public SectionClass section;
        public MaterialClass material;
        public Line axis;
        public NodeClass startNode;
        public NodeClass endNode;
        


        public BeamClass()
        { }
        
        public BeamClass(string _Name, Line _line, SectionClass _section, MaterialClass _material) 
        {
            Name = _Name;
            axis= _line;
            section = _section;
            material = _material;
        }
    }
}
