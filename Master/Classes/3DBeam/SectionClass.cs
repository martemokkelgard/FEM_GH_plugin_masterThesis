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
        public double w;
        public double h;
        public double tw;
        public double tf;
        public double r;
        public double CSA;
        public double Iy;
        public double Iz;
        public double Iyz;
        public double J;
        public double a;
        public double b;
        public double Asteg;



        public SectionClass()
        {
            //empty
        }

        public SectionClass(string _Name, double _h, double _w, double _tf, double _tw)
        {

            Name = _Name;
            w = _w;
            h = _h;
            tw = _tw;
            tf = _tf;
            CSA = w * h - ((w - 2 * tw) * (h - 2 * tf));
            Iy = (1 / (12.00) * w * Math.Pow(h, 3)) - (1 / (12.00) * (w - 2 * tw) * Math.Pow((h - 2 * tf), 3));
            Iz = (1 / (12.00) * h * Math.Pow(w, 3)) - (1 / (12.00) * (h - 2 * tf) * Math.Pow((w - 2 * tw), 3));
            Asteg = (h - 2 * tf) * 2;


            if (h > w)
            {
                a = h / 2;
                b = w / 2;

                J = a * Math.Pow(b, 3) * (16.00 / 3.00 - 3.36 * b / a * (1 - Math.Pow(b, 4) / (12.00 * Math.Pow(a, 4)))); //polar moment of interia
            }

            else
            {
                b = h / 2;
                a = w / 2;
                J = a * Math.Pow(b, 3) * (16.00 / 3.00 - 3.36 * b / a * (1 - Math.Pow(b, 4) / (12.00 * Math.Pow(a, 4))));
            }



        }

        public SectionClass(string _Name, double _r)
        {
            Name = _Name;
            r = _r;
            CSA = Math.PI * Math.Pow(r, 2);
            Iy = (Math.PI * Math.Pow(r, 4)) / 4.0;
            Iz = Iy;
            J = (Math.PI * Math.Pow(r, 4)) / 2.0;

        }
    }
}


