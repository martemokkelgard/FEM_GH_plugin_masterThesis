using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BcClass
    {

        public Point3d Coordinate;
        public bool ux;
        public bool uy;
        public bool uz;
        public bool rx;
        public bool ry;
        public bool rz;

        public BcClass()
        { }



        public BcClass(Point3d _Coordinate,bool _ux, bool _uy, bool _uz, bool _rx, bool _ry, bool _rz)
        {
            
            Coordinate = _Coordinate;
            ux = _ux;
            uy = _uy;
            uz = _uz;
            rx = _rx;
            ry = _ry;
            rz = _rz;
        }

    }
}
