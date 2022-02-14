using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class CreatePointLoad : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Load class.
        /// </summary>
        public CreatePointLoad()
          : base("PointLoad", "Nickname",
              "Description",
              "Løve", "2DTruss")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply load(s)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Load", "L", "Load magnitude [Newtons].Give either one load to be applied to all inputted points, or different loads for each inputted loads", GH_ParamAccess.list);
            pManager.AddNumberParameter("angle (xz)", "a", "Angle [degrees] for load in xz plane, default: 90deg", GH_ParamAccess.list, 90);  
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("PointLoad", "PL", "List of PointLoadsClass object", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            
            //input
            List<Point3d> pointList = new List<Point3d>(); //list of points where loads are applied
            List<double> AmpList = new List<double>();  //list of amplitudes of the loads
            List<double> anglexzList = new List<double>();  //list of angles in xz plane of the loads (2D)
            List<Vector3d> Vecs = new List<Vector3d>();

            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetDataList(1, AmpList)) return;
            DA.GetDataList(2, anglexzList);

            //code
            
            double load = 0;
            double vecX = 0;
            double vecY = 0;
            double vecZ = 0;

            if (AmpList.Count == 1 && anglexzList.Count == 1)  //similar load amplitude and angle for all points
            {
                
                load = -1 * AmpList[0]; //neg for z-dir
                vecX = Math.Round(load * Math.Cos(anglexzList[0] * Math.PI / 180) * Math.Cos(anglexzList[0] * Math.PI / 180), 2);
                vecY = Math.Round(load * Math.Cos(anglexzList[0] * Math.PI / 180) * Math.Sin(anglexzList[0] * Math.PI / 180), 2);
                vecZ = Math.Round(load * Math.Sin(anglexzList[0] * Math.PI / 180), 2);

                Vector3d xyzVec = new Vector3d(vecX,vecY,vecZ);
                

                for (int i = 0; i < pointList.Count; i++)
                {
                    Vecs.Add(xyzVec);
                }

            }

            else  //different load amplitude and angle for the points
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    if (AmpList.Count < i) //loadlist is shorter than pointlist
                    {
                        Vector3d xyzVec = new Vector3d(vecX, vecY, vecZ);
                        Vecs.Add(xyzVec);
                    }
                    else
                    {
                        
                        load = -1 * AmpList[0]; //neg for z-dir
                        vecX = Math.Round(load * Math.Cos(anglexzList[i]* Math.PI / 180) * Math.Cos(anglexzList[i]*Math.PI / 180), 2);
                        vecY = Math.Round(load * Math.Cos(anglexzList[i] * Math.PI / 180) * Math.Sin(anglexzList[i] * Math.PI / 180), 2);
                        vecZ = Math.Round(load * Math.Sin(anglexzList[i] * Math.PI / 180), 2);

                        Vector3d xyzVec = new Vector3d(vecX, vecY, vecZ);
                        Vecs.Add(xyzVec);
                    }

                    
                }
            
            
            }
            List<LoadClass> loads = new List<LoadClass>();
            for (int i = 0; i < pointList.Count; i++)
            {
                loads.Add(new LoadClass(pointList[i], Vecs[i]));
            }



            //output
            DA.SetDataList(0, loads);

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("ED282DDC-4338-474F-BBB1-0B01FD70D9F5"); }
        }
    }
}