using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BCClass2D
    {

        public Point3d Coordinate;
        public bool ux;
        public bool uz;
        public bool rx;

        public BCClass2D()
        { }



        public BCClass2D(Point3d _Coordinate, bool _ux, bool _uz, bool _rx)
        {

            Coordinate = _Coordinate;
            ux = _ux;
            uz = _uz;
            rx = _rx;
        }

    }
}