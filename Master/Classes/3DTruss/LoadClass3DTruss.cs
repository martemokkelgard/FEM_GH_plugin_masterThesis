using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class LoadClass3DTruss
    {
        public Point3d coordinate;
        public Vector3d LoadVec;

        public LoadClass3DTruss()
        { }


        public LoadClass3DTruss(Point3d _coordinate, Vector3d _LoadVec)
        {
            coordinate = _coordinate;
            LoadVec = _LoadVec;

        }



    }
}
