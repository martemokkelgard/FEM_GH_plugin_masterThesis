using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class BC3DBeam : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the BoundaryConditions class.
        /// </summary>
        public BC3DBeam()
          : base("BC3DBeam", "Nickname",
              "Description",
              "Panda", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points","P", "Points for BCs", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Ux", "Ux", "Ux: false = free", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Uy", "Uy", "Uy: false = free", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Uz", "Uz", "Uz: false = free", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Rx", "Rx", "Rx: false = free", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Ry", "Ry", "Ry: false = free", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Rz", "Rz", "Rz: false = free", GH_ParamAccess.item, false);
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
            bool y = false;
            bool z = false;
            bool rx = false;
            bool ry = false;
            bool rz = false;


            if (!DA.GetDataList(0, pt)) return;
            if (!DA.GetData(1, ref x)) return;
            if (!DA.GetData(2, ref y)) return;
            if (!DA.GetData(3, ref z)) return;
            if (!DA.GetData(4, ref rx)) return;
            if (!DA.GetData(5, ref ry)) return;
            if (!DA.GetData(6, ref rz)) return;

            List<BcClass> BCs = new List<BcClass>();
          

            foreach (Point3d p in pt)
            {
                BCs.Add(new BcClass(p, x, y, z, rx, ry, rz ));
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
            get { return new Guid("af95171f-3b5c-4e23-866d-58642e5021b2"); }
        }
    }
}