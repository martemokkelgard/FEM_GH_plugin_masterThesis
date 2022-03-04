using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System;
using System.Collections.Generic;
using System.Linq;
     


namespace Master
{
    public class MasterComponent : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public MasterComponent()
          : base("Master", "Nickname",
            "2DTruss FEM",
            "Løve", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
           //pManager.AddLineParameter("Lines", "LNS", "Geometry (lines)", GH_ParamAccess.list);
            //pManager.AddTextParameter("Boundary Conditions", "BDC", "Boundary conditions ((x,z:1,1 where 1 = free and 0 = fixed))",GH_ParamAccess.list);
            //pManager.AddNumberParameter("Crossection Area", "A", "Area of cross section, initial value 10e3", GH_ParamAccess.item, 10000);
            //pManager.AddNumberParameter("E modulus of material")
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// You can add image files to your project resources and access them like this:
        /// return Resources.IconForThisComponent;
        /// </summary>
        protected override System.Drawing.Bitmap Icon => null;

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid => new Guid("38C36126-5A67-49E6-B022-5CD0F28B2DC7");
    }
}