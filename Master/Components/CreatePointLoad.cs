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
            pManager.AddVectorParameter("Load", "L", "Load vector with amplitude in [N].Give either one load to be applied to all inputted points, or different loads for each inputted point", GH_ParamAccess.item);
            
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
            Vector3d Vecs = new Vector3d();

            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetData(1, ref Vecs)) return;
            

            //code
            
            
            List<LoadClass> loads = new List<LoadClass>();
            for (int i = 0; i < pointList.Count; i++)
            {
                loads.Add(new LoadClass(pointList[i], Vecs));
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