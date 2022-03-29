using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class CreateBar : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Beam class.
        /// </summary>
        public CreateBar()
          : base("CreateBar", "Nickname",
              "Description",
              "Løve", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Name", "N", "Name of bar", GH_ParamAccess.item);
            pManager.AddCurveParameter("lines", "LNS", "Geometry (lines with dim. in [m])", GH_ParamAccess.list);
            pManager.AddGenericParameter("material", "M", "MaterialClass", GH_ParamAccess.item);
            pManager.AddGenericParameter("section", "S", "Section of the bar", GH_ParamAccess.item);
            
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("bars", "B", "List of barClass objects", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //input
            string name = "Bar";
            List<Curve> lines = new List<Curve>();
            MaterialClass mat = new MaterialClass();
            SectionClass sec = new SectionClass();
            DA.GetData(0, ref name);
            DA.GetDataList(1, lines);
            DA.GetData(2, ref mat);
            DA.GetData(3, ref sec);
            List<BarClass> bars = new List<BarClass>();

            
            //code
            foreach (Curve l in lines)  //making barsClass objects of lines
            { 
                bars.Add(new BarClass("trussBar", l, sec, mat));
            }

            for (int i = 0; i < bars.Count; i++)   //giving id to beamClass objects
            {
                bars[i].Id = i;
            }

           
            List<Point3d> pts = new List<Point3d>();  //making pointList from lines

            for (int i = 0; i < lines.Count; i++)
            {
                Curve L1 = lines[i];

                if (!pts.Contains(L1.PointAtStart))
                {
                    pts.Add(L1.PointAtStart);
                }
                if (!pts.Contains(L1.PointAtEnd))
                {
                    pts.Add(L1.PointAtEnd);
                }

            }

            List<NodeClass> nodes = new List<NodeClass>();  //making nodeClass objects of points
            foreach (Point3d p in pts)
            {
                nodes.Add(new NodeClass(p));
            }

            
            for (int i = 0; i < nodes.Count; i++)   //giving id to nodeClass objects
            {
                nodes[i].Id = i;

                foreach (BarClass l in bars)       //finding the nodeClass object that is start/end node of barClass objects
                {


                    if (l.axis.PointAtStart == nodes[i].pt)
                    {
                        l.startNode = nodes[i];
                    }
                    if (l.axis.PointAtEnd == nodes[i].pt)
                    {
                        l.endNode = nodes[i];
                    }


                }

            }



            
            //output
            DA.SetDataList(0, bars);
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
            get { return new Guid("A5FAD6C5-F8A4-44BC-A8D3-FDF88A0CB043"); }
        }
    }
}