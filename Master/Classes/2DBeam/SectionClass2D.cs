using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class SectionClass2D
    {
        public string Name;
        public int id;
        public double w;
        public double h;
        public double tw;
        public double tf;
        public double CSA;
        public double I;


        public SectionClass2D()
        {
            //empty
        }

        public SectionClass2D(string _Name, double _h, double _w, double _tf, double _tw)
        {

            Name = _Name;
            w = _w;
            h = _h;
            tw = _tw;
            tf = _tf;
            CSA = w * h - ((w - 2 * tw) * (h - 2 * tf));
            I = (1 / (12.0000) * w * Math.Pow(h, 3)) - (1 / (12.0000) * (w - 2 * tw) * Math.Pow((h - 2 * tf), 3));

        }
    }
}