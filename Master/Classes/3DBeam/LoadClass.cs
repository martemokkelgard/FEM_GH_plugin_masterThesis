using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master
{
    internal class LoadClass
    {
        public Point3d coordinate;
        public Vector3d LoadVec;
        public Curve line;
        public Point3d stPt;
        public Point3d enPt;
        public bool Id;
        public bool Id1;
        public List<double> mVal;

        public LoadClass()
        { }


        public LoadClass(Point3d _coordinate, Vector3d _LoadVec)
        {
           
            coordinate = _coordinate;
            LoadVec = _LoadVec;
            Id = true;
            Id1 = true; 
        }

        public LoadClass(Curve _line, Vector3d _LoadVec)
        {
            line = _line;
            LoadVec = _LoadVec/2;
            stPt = line.PointAtStart;
            enPt = line.PointAtEnd;
            Id = false;
        }

        public LoadClass(Point3d _coordinate, List<double> _momentValues)
        {
            coordinate = _coordinate;
            Id = true;
            mVal = _momentValues;
            Id1 = false;
        }

    }
}
