using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BcClass2DTruss
    {

        public Point3d Coordinate;
        public bool ux;
        public bool uz;

        public BcClass2DTruss()
        { }



        public BcClass2DTruss(Point3d _Coordinate, bool _ux, bool _uz)
        {

            Coordinate = _Coordinate;
            ux = _ux;
            uz = _uz;
        }

    }
}