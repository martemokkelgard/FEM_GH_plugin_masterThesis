using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class Section : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Sectioncs class.
        /// </summary>
        public Section()
          : base("Sectioncs", "Nickname",
              "Description",
              "Løve", "3DTruss")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "N", "Name of the section (bxh [mm]), default: 100x100mm", GH_ParamAccess.item,"100x100");
            pManager.AddNumberParameter("Height", "h", "Height of the cross-section", GH_ParamAccess.item, 100);
            pManager.AddNumberParameter("Width", "w", "Width of the cross-section", GH_ParamAccess.item, 100);
            pManager.AddNumberParameter("ThicknessF", "tf", "Thickness of the flanges of the cross-section", GH_ParamAccess.item, 5);
            pManager.AddNumberParameter("ThicknessW", "tw", "Thickness of the web of the cross-section", GH_ParamAccess.item, 5);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("section","S","SectionClass object",GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            //input
            string Name = "100x100";
            double h = 100;
            double w = 100;
            double tf = 5;
            double tw = 5;
            DA.GetData(0, ref Name);
            DA.GetData(1,ref h);
            DA.GetData(2,ref w);
            DA.GetData(3,ref tf);
            DA.GetData(4,ref tw);
               
            //code
            SectionClass sec = new SectionClass(Name, h, w, tf, tw);
           

            //output
            DA.SetData(0, sec);
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
            get { return new Guid("D6B5319A-CE46-464C-ABFF-EC3CCC6DFB5B"); }
        }
    }
}