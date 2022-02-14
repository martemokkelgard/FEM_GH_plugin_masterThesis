using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class SectionClass
    {
        public string Name;
        public int id;
        public double width;
        public double height;
        public double CSA;


        public SectionClass()
        {
            //empty
        }
        
        public SectionClass(string _Name)
        {
            Name = _Name;

            string[] dimensions = Name.Split('x');
            width = Convert.ToDouble(dimensions[0]);
            height = Convert.ToDouble(dimensions[1]);
            CSA = width * height;
                   
        }
    }
}

