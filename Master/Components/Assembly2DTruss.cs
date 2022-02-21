using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;

//hei

namespace Master.Components
{
    public class Assembly2DTruss : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Assembly2DTruss class.
        /// </summary>
        public Assembly2DTruss()
          : base("Assembly2DTruss", "Nickname",
              "Description",
              "Løve", "2DTruss")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Bar", "B", "BarClass object", GH_ParamAccess.list);
            pManager.AddGenericParameter("Boundary Conditions", "BC", "BcClass object", GH_ParamAccess.list);
            pManager.AddGenericParameter("Point load", "PL", "LoadClass object", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Defomation", "def", "Deformations of the bar in [mm]", GH_ParamAccess.list);
            pManager.AddNumberParameter("ReactionForces", "R", "Reaction Forces", GH_ParamAccess.list);
            pManager.AddPointParameter("DisplacementOfNodes", "displ", "Deformed geometry", GH_ParamAccess.list);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            
            
            //input
            
            
            List<Point3d> BcPts = new List<Point3d>();  //hei
            List<bool> BcValue = new List<bool>();
            List<Point3d> LoadPts = new List<Point3d>();
            List<Vector3d> LoadVec = new List<Vector3d>();

            MaterialClass mat = new MaterialClass();
            SectionClass sec = new SectionClass();
            List<BarClass> bars = new List<BarClass>();
            List<LoadClass> lc = new List<LoadClass>();
            List<BcClass> bcc = new List<BcClass>();

            
            //int node1 = bars[1].endNode.Id;

            DA.GetDataList(0, bars);
            DA.GetDataList(1,  bcc);
            DA.GetDataList(2, lc);

            /*
            double E = bars[i].material.youngsModolus;
            double A = bars.section.CSA;
            lines = bc.lines;
            BcPts = bcc.Coordinates;
            BcValue = bcc.xyzVec;
            LoadPts = lc.coordinate;
            LoadAmplitude = lc.xyzVec;
            */

            foreach (var load in lc)
            {
                LoadPts.Add(load.coordinate);
                LoadVec.Add(load.LoadVec);
            }


            //code

            List<Point3d> pts = new List<Point3d>();
            foreach (BarClass b in bars)
            {
                
                if (!pts.Contains(b.startNode.pt))
                {
                    pts.Add(b.startNode.pt);
                }
                if (!pts.Contains(b.endNode.pt))
                {
                    pts.Add(b.endNode.pt);
                }
            }
            
            List<int> BCList = CreateBCList(bcc, pts);
            List<double> LoadList = CreateLoadList(LoadPts, LoadVec, pts);
            List<int> BC = new List<int>();

            

           
            int dofs_red = BC.Sum();

            //Vector<double> R_red = SparseVector.OfEnumerable(new double[dofs_red]);
            Matrix<double> K_red = SparseMatrix.OfArray(new double[dofs_red, dofs_red]);
            //List<double> ptsK = new List<double>();
            
            Matrix<double> k_eG = DenseMatrix.OfArray(new double[2, 2]);
            Matrix<double> k_tot = DenseMatrix.OfArray(new double[pts.Count*2, pts.Count*2]);
            Vector<double> forces = SparseVector.OfEnumerable(new double[4]);
            



            CreateGlobalStiffnesMatrix(bars, pts, out k_tot, out k_eG);

            

            
            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(BCList, k_tot, out K_red);

            var invK = K_red.Inverse();
            
            var def = invK.Multiply(R);

            var displNodes = new List<Point3d>();


            CreateForces(bars, pts, def, out forces);


            for (int i =0; i< pts.Count; i++)
            {

                
                displNodes.Add(new Point3d(pts[i].X + def[2 * i]/1000, pts[i].Y , pts[i].Z + def[2 * i + 1]/1000));
            }


            DA.SetDataList(0, def);
            DA.SetDataList(0, forces);
            DA.SetDataList(2, displNodes);
            //output

        }

        private List<int> CreateBCList( List<BcClass> _BcValue, List<Point3d> _Pts)  //making a list with indexes of fixed BC
        {
            
            //List<int> BCs = new List<int>();
            List<int> BCsIndex = new List<int>();

            
            for (int i = 0; i < _Pts.Count; i++)
            {
                Point3d node = _Pts[i];
                foreach (var b in _BcValue)
                {
                    if (b.Coordinate.DistanceTo(node) < 0.000001)
                    {
                        if(b.ux)
                            BCsIndex.Add(2*i);

                        if(b.uz)
                            BCsIndex.Add(2*i+1);
                    }
                }

            }

            return BCsIndex;
        }

        private List<double> CreateLoadList(List<Point3d> _LoadPts, List<Vector3d> _LoadValue, List<Point3d> _Pts)
        {
            
            List<double> LoadValue = new List<double>();


            for (int i = 0; i < _Pts.Count; i++)
            {
                //this line was not so good
                Point3d node = _Pts[i];


                for (int j = 0; j < _LoadPts.Count; j++)
                {
                    Point3d loadLocation = _LoadPts[j];
                    if (loadLocation.DistanceTo(node) < 0.00001)
                    {
                        LoadValue.Add(_LoadValue[j][0]);
                        LoadValue.Add(_LoadValue[j][2]);
                    }
                    else
                    {
                        LoadValue.Add(0);
                        LoadValue.Add(0);
                    }

                }

            }

            return LoadValue;
        }

        private static void CreateForces(List<BarClass> bars, List<Point3d> points, Vector<double> _def, out Vector<double> forces)
        {
            Matrix<double> k_eG = DenseMatrix.OfArray(new double[2, 2]);
            Vector<double> S = SparseVector.OfEnumerable(new double[4]);
            Vector<double> v = SparseVector.OfEnumerable(new double[4]);

            foreach (BarClass b in bars)
            {
                Line currentLine = b.axis;
                double mat = (b.section.CSA * b.material.youngsModolus) / (currentLine.Length * 1000);
                Point3d p1 = currentLine.From;
                Point3d p2 = currentLine.To;

                double angle = Math.Atan2(p2.Z - p1.Z, p2.X - p1.X);  //returns angle in rad
                double c = Math.Cos(angle); //input angle should be in rad
                double s = Math.Sin(angle);

                Matrix<double> T = SparseMatrix.OfArray(new double[,]
                                    {
                        { c, s, 0 ,0},
                        { -s, c ,0, 0},
                        { 0, 0, c, -s},
                        { 0 ,0 ,s, c}
                                    });

                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                                    {
                        { 1, 0, -1 ,0},
                        { 0, 0 ,0, 0},
                        { -1, 0, 1, 0},
                        { 0 ,0 ,0, 0}
                                    });

                Matrix<double> Tt = T.Transpose(); //transpose
                k_eG = Tt.Multiply(ke);
                k_eG = k_eG.Multiply(T);


                k_eG = mat * k_eG;  //global element stivehetsmatrise
                ke = mat * ke;

                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                v[0] = _def[node1 * 2];
                v[1] = _def[node1 * 2 +1];
                v[2] = _def[node2 * 2];
                v[3] = _def[node2 * 2 +1];

                S = k_eG.Multiply(v);


            }

            forces = S;
        }

        private static void CreateGlobalStiffnesMatrix(List<BarClass> bars, List<Point3d> points, out Matrix<double> k_tot ,out Matrix<double> k_eG)

        {

            int dofs = points.Count * 2;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[dofs, dofs]);
            Matrix<double> K_eG = DenseMatrix.OfArray(new double[2, 2]);

            foreach (BarClass b in bars)
            {
                            

                Line currentLine = b.axis;
                double mat = (b.section.CSA * b.material.youngsModolus) / (currentLine.Length*1000);
                Point3d p1 = currentLine.From;
                Point3d p2 = currentLine.To;

                double angle = Math.Atan2(p2.Z - p1.Z, p2.X - p1.X);  //returns angle in rad
                double c = Math.Cos(angle); //input angle should be in rad
                double s = Math.Sin(angle);

                /*
                double[,] K_e = new double[,]
                {   { mat*c*c, mat*c*s, -mat*c*c ,-mat*c*s},
                        { mat*c*s, mat*s*s ,-mat*c*s, -mat*s*s},
                        { -mat*c*c, -mat*c*s, mat*c*c, mat*c*s},
                        { -mat*c*s ,-mat*s*s ,mat*c*s, mat*s*s} };
                */
                Matrix<double> T = SparseMatrix.OfArray(new double[,]
                                {   
                        { c, s, 0 ,0},
                        { -s, c ,0, 0},
                        { 0, 0, c, -s},
                        { 0 ,0 ,s, c} 
                                });

                Matrix<double> ke = SparseMatrix.OfArray(new double[,]
                                {
                        { 1, 0, -1 ,0},
                        { 0, 0 ,0, 0},
                        { -1, 0, 1, 0},
                        { 0 ,0 ,0, 0}
                                });



                Matrix<double> Tt = T.Transpose(); //transpose
                k_eG = Tt.Multiply(ke);
                K_eG = k_eG.Multiply(T);  
                

                K_eG = mat * K_eG;  //global element stivehetsmatrise
                ke = mat * ke;      //lokal element stivhetsmatrise

                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                //upper left corner of k-matrix
                K_tot[node1 * 2, node1 * 2] += K_eG[0,0];
                K_tot[node1 * 2, node1 * 2 + 1] += K_eG[0, 1];
                K_tot[node1 * 2 + 1, node1 * 2] += K_eG[1, 0];
                K_tot[node1 * 2 + 1, node1 * 2 + 1] += K_eG[1, 1];
                
                //upper right corner of k-matrix
                K_tot[node1 * 2, node2 * 2] += K_eG[0, 2];
                K_tot[node1 * 2, node2 * 2 + 1] += K_eG[0, 3];
                K_tot[node1 * 2 + 1, node2 * 2] += K_eG[1, 2];
                K_tot[node1 * 2 + 1, node2 * 2 + 1] += K_eG[1, 3];
                
                //lower left corner of k-matrix
                K_tot[node2 * 2, node1 * 2] += K_eG[2, 0];
                K_tot[node2 * 2, node1 * 2 + 1] += K_eG[2, 1];
                K_tot[node2 * 2 + 1, node1 * 2] += K_eG[3, 0];
                K_tot[node2 * 2 + 1, node1 * 2 + 1] += K_eG[3, 1];
                
                 //lower right corner of k-matrix
                K_tot[node2 * 2, node2 * 2] += K_eG[2, 2];
                K_tot[node2 * 2, node2 * 2 + 1] += K_eG[2, 3];
                K_tot[node2 * 2 + 1, node2 * 2] += K_eG[3, 2];
                K_tot[node2 * 2 + 1, node2 * 2 + 1] += K_eG[3, 3];
            

            }

            k_eG = K_eG;
            k_tot = K_tot;


        }

        private static void CreateReducedStiffnesMatrix(List<int> _BC, Matrix<double> K_tot, out Matrix<double> K_red)
        {

            
            for (int i = 0; i < _BC.Count; i++)
            {
                int sizeOfMatrix = K_tot.ColumnCount;
                int blockedIndex = _BC[i];

            

                for (int r = 0; r < sizeOfMatrix; r++)
                {
                    K_tot[blockedIndex, r] = 0;
                    K_tot[r, blockedIndex] = 0;

                }


                K_tot[blockedIndex, blockedIndex]=1;


            }

            K_red = K_tot;
            
          

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
            get { return new Guid("DDB74701-399A-4352-B4DC-1F175E566460"); }
        }
    }
}