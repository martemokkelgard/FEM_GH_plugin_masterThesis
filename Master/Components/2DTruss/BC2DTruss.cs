using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components._2DTruss
{
    public class BC2DTruss : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the BC2DTruss class.
        /// </summary>
        public BC2DTruss()
         : base("BC2DTruss", "Nickname",
             "Description",
             "Panda", "2DTruss")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points for BCs", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Ux", "Ux", "Ux", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Uz", "Uz", "Uz", GH_ParamAccess.item, false);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Boundary Conditions", "BCs", "BcClass object", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            List<Point3d> pt = new List<Point3d>();
            bool x = false;
            bool z = false;

            if (!DA.GetDataList(0, pt)) return;
            if (!DA.GetData(1, ref x)) return;
            if (!DA.GetData(2, ref z)) return;

            List<BcClass2DTruss> BCs = new List<BcClass2DTruss>();


            foreach (Point3d p in pt)
            {
                BCs.Add(new BcClass2DTruss(p, x, z));
            }

            //output
            DA.SetDataList(0, BCs);
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
            get { return new Guid("F901549A-C8B9-46D5-B96D-B8C3C5C35D89"); }
        }
    }
}