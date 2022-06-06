using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class CreateMomentForce : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateMomentForce class.
        /// </summary>
        public CreateMomentForce()
          : base("CreateMomentForce", "Nickname",
              "Description",
              "Panda", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply load(s)", GH_ParamAccess.list);
            pManager.AddNumberParameter("MomentValue", "MV", "Moment Load vector with amplitude in [Nmm].Give it like list", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("MomentLoad", "ML", "List of MomentLoadsClass object", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            List<Point3d> pointList = new List<Point3d>(); //list of points where loads are applied
            List<double> mVal = new List<double>();


            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetDataList(1, mVal)) return;



            //code


            List<LoadClass> loads = new List<LoadClass>();

            for (int i = 0; i < pointList.Count; i++)
            {
                loads.Add(new LoadClass(pointList[i], mVal));
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
            get { return new Guid("83489457-ad76-4d3e-b757-22f1255d3aa8"); }
        }
    }
}