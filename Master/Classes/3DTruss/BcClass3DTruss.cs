using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class BcClass3DTruss
    {

        public Point3d Coordinate;
        public bool ux;
        public bool uy;
        public bool uz;

        public BcClass3DTruss()
        { }



        public BcClass3DTruss(Point3d _Coordinate, bool _ux, bool _uy, bool _uz)
        {

            Coordinate = _Coordinate;
            ux = _ux;
            uy = _uy;
            uz = _uz;
        }

    }
}
