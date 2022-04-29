using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Master.Components
{
    
    public class CreateBeamNonlinear : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Beam class.
        /// </summary>
        public CreateBeamNonlinear()
          : base("CreateBeamNonlinear", "Nickname",
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
            pManager.AddNumberParameter("Degree", "Deg", "Degree of nonlinearity", GH_ParamAccess.item, 1);

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
            /*
            //input
            string name = "Beam";
            List<Curve> lines = new List<Curve>();
            MaterialClass mat = new MaterialClass();
            SectionClass sec = new SectionClass();
            double deg = 1;
            DA.GetData(0, ref name);
            DA.GetDataList(1, lines);
            DA.GetData(2, ref mat);
            DA.GetData(3, ref sec);
            DA.GetData(4, ref deg);
            List<BeamClassNonLin> beams = new List<BeamClassNonLin>();

            
            lines[0].PointAtNormalizedLength(0);
            lines[0].PointAtNormalizedLength(0.5);
            lines[0].PointAtNormalizedLength(1);
            double minpdomain = lines[0].Domain[0];
            double maxpdomain = lines[0].Domain[1];
            double middomain = minpdomain + (maxpdomain - minpdomain) / 2;
            lines[0].PointAt(minpdomain);
            lines[0].PointAt(middomain);
            lines[0].PointAt(maxpdomain);
            


            //code
            foreach (Curve l in lines)  //making barsClass objects of lines
            {
                beams.Add(new BeamClassNonLin("Beam", l, sec, mat));
            }

            for (int i = 0; i < beams.Count; i++)   //giving id to beamClass objects
            {
                beams[i].Id = i;
            }

            List<double> param = new List<double>();
            for (int i = 0; i < deg + 1; i++)
            {
                param.Add(i / deg);
            }



            List<Point3d> pts = new List<Point3d>();  //making pointList from lines

            for (int i = 0; i < lines.Count; i++)
            {
                Curve L1 = lines[i];
                foreach (double p in param)
                {
                    pts.Add(L1.PointAtNormalizedLength(p));
                }

            }

            for (int i = 0; i < pts.Count; i++)
            {
                for (int j = 0; j < pts.Count; j++)
                {
                    double tol = 0.00001;
                    if (pts[i].DistanceTo(pts[j]) < tol && i != j)
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


                foreach (BeamClassNonLin l in beams)       //finding the nodeClass object that is start/end node of barClass objects
                {

                    double tol = 0.001;

                    foreach (double p in param)
                    {
                        if (l.axis.PointAtNormalizedLength(p).DistanceTo(nodes[i].pt) < tol)
                        {
                             
                            int u = ((int)(p * deg));
                            l.nodes[u] = nodes[i];
                        }
                    }



                }

            }




            //output
            DA.SetDataList(0, beams);
            */
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
            get { return new Guid("636494A3-78AA-42DE-BBF0-7FC418618C2D"); }
        }
    }
}