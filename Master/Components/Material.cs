using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class Material : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Material class.
        /// </summary>
        public Material()
          : base("Material", "Nickname",
              "Description",
              "Panda", "3DTruss")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "N", "name of the materiaL, default: S355", GH_ParamAccess.item, "s355");
            pManager.AddNumberParameter("Density", "D", "density of the material, default: 7.8e-6 kg/mm^3", GH_ParamAccess.item, 0.0000078 );
            pManager.AddNumberParameter("YoungsModulus", "E", "youngsModulus of the material (in [N/mm^2], default 210000 N/mm^2", GH_ParamAccess.item,210000);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("material", "M", "materialClass object", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            string Name = "C14";
            double d = 100;
            double E = 100000;

            DA.GetData(0, ref Name);
            DA.GetData(1, ref d);
            DA.GetData(2, ref E);

            //code
            MaterialClass mat = new MaterialClass(Name, d, E);

            //output
            DA.SetData(0, mat);

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
            get { return new Guid("0e3b608a-f51e-44a6-9bd0-6671f57735ae"); }
        }
    }
}