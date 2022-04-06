using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class CreateLineLoad : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the CreateLineLoad class.
        /// </summary>
        public CreateLineLoad()
          : base("CreateLineLoad", "Nickname",
              "Description",
              "Panda", "2DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddLineParameter("Line", "L", "Lines to add distributed load to", GH_ParamAccess.list);
            pManager.AddVectorParameter("LoadVector", "LV", "Load vector with amplitude in [N].Give either one load to be applied to all lines", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("LineLoad", "PL", "List of PointLoadsClass object", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            Vector3d Vecs = new Vector3d();
            List<Line> lines = new List<Line>();

            if (!DA.GetDataList(0, lines)) return;
            if (!DA.GetData(1, ref Vecs)) return;

            //code

            List<LoadClass> loads = new List<LoadClass>();

            for (int i = 0; i < lines.Count; i++)
            {
                loads.Add(new LoadClass(lines[i], Vecs));

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
            get { return new Guid("1e7fd7a0-0fc1-47a1-8dd8-af4c50507d5f"); }
        }
    }
}