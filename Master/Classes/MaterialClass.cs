using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class MaterialClass
    {
        public string Name;
        public int Id;
        public double youngsModolus;
        public double density;
        public double G;

        public MaterialClass()
        { }

        public MaterialClass(string _Name)
        {
            Name = _Name;
        }


        public MaterialClass(string _Name, double _density, double _youngsModolus)
        {
            Name = _Name;
            density = _density;
            youngsModolus = _youngsModolus;
            G = 1;
        }


    }
}
