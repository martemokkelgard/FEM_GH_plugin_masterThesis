using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    public class CreateBarBilinear : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Beam class.
        /// </summary>
        public CreateBarBilinear()
          : base("CreateBarBilinear", "Nickname",
              "Description",
              "Panda", "3DBeam")
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

            /*
            lines[0].PointAtNormalizedLength(0);
            lines[0].PointAtNormalizedLength(0.5);
            lines[0].PointAtNormalizedLength(1);

            double minpdomain = lines[0].Domain[0];
            double maxpdomain = lines[0].Domain[1];
            double middomain = minpdomain + (maxpdomain - minpdomain) / 2;
            lines[0].PointAt(minpdomain);
            lines[0].PointAt(middomain);
            lines[0].PointAt(maxpdomain);
            */


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
                
                double a = 0.5;                                             // * parameter for midnode 

                pts.Add(new Point3d(Math.Round(L1.PointAtNormalizedLength(0).X,6), Math.Round(L1.PointAtNormalizedLength(0).Y, 6), Math.Round(L1.PointAtNormalizedLength(0).Z, 6)));

                pts.Add(new Point3d(Math.Round(L1.PointAtNormalizedLength(a).X, 6), Math.Round(L1.PointAtNormalizedLength(a).Y, 6), Math.Round(L1.PointAtNormalizedLength(a).Z, 6)));

                pts.Add(new Point3d(Math.Round(L1.PointAtNormalizedLength(1).X, 6), Math.Round(L1.PointAtNormalizedLength(1).Y, 6), Math.Round(L1.PointAtNormalizedLength(1).Z, 6)));

            }

            for (int i = 0;i < pts.Count;i++)
            {
                for (int j = 0;j < pts.Count;j++)
                {
                    double tol = 0.00001;
                    if (pts[i].DistanceTo(pts[j]) < tol && i!=j)
                    {
                        pts.Remove(pts[j]);
                    }
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
                    
                    double a = 0.5;                                             // * parameter for midnode 
                    double tol = 0.001;

                    if (l.axis.PointAtStart.DistanceTo(nodes[i].pt)<tol)
                    {
                        l.startNode = nodes[i];
                    }


                    if (l.axis.PointAtNormalizedLength(a).DistanceTo(nodes[i].pt)< tol)
                    {
                        l.midNode = nodes[i];
                    }

                    

                    if (l.axis.PointAtEnd.DistanceTo(nodes[i].pt)< tol)
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
            get { return new Guid("fe5563f1-5582-4065-a219-fb8ffa5aaec0"); }
        }
    }
}