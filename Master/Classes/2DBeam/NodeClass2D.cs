using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class NodeClass2D
    {
        public int Id;
        public Line axis;
        public Point3d pt;

        public NodeClass2D()
        {

        }

        public NodeClass2D(Point3d _pt) //use barclass.axis
        {
            pt = _pt;



        }
    }


}