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
            
            
            List<Point3d> BcPts = new List<Point3d>();
            List<bool> BcValue = new List<bool>();
            List<Point3d> LoadPts = new List<Point3d>();
            List<Vector3d> LoadVec = new List<Vector3d>();

            MaterialClass mat = new MaterialClass();
            SectionClass sec = new SectionClass();
            List<BarClass> bars = new List<BarClass>();
            List<LoadClass> lc = new List<LoadClass>();
            List<BcClass> bcc = new List<BcClass>();

            //int node1 = bars[i].endNode.Id

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

            Vector<double> R_red = SparseVector.OfEnumerable(new double[dofs_red]);
            Matrix<double> K_red = SparseMatrix.OfArray(new double[dofs_red, dofs_red]);
            List<double> ptsK = new List<double>();
            
            Matrix<double> k_eG = SparseMatrix.OfArray(new double[2, 2]);
            Matrix<double> k_tot = SparseMatrix.OfArray(new double[pts.Count*2, pts.Count*2]);
            

            CreateGlobalStiffnesMatrix(bars, pts, out k_tot, out k_eG);


            
            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(pts,BCList, k_tot, R, out K_red, out R_red, out ptsK);

            var invK = K_red.Inverse();
            
            var def = invK.Multiply(R_red);

            var displNodes = new List<Point3d>();

            foreach (BarClass b in bars)
            {
                //var S = k_eG * def[];
            }

            

            
            for(int i =0; i< ptsK.Count/2; i++)
            {

                
                displNodes.Add(new Point3d(ptsK[i] + def[2 * i]/1000, 0, ptsK[i] + def[2 * i + 1]/1000));
            }


            DA.SetDataList(0, def);
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


        private static void CreateGlobalStiffnesMatrix(List<BarClass> bars, List<Point3d> points, out Matrix<double> k_tot ,out Matrix<double> k_eG)

        {

            int dofs = points.Count * 2;
            Matrix<double> K_tot = SparseMatrix.OfArray(new double[dofs, dofs]);
            Matrix<double> K_eG = SparseMatrix.OfArray(new double[2, 2]);

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

                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                                {
                        { 1, 0, -1 ,0},
                        { 0, 0 ,0, 0},
                        { -1, 0, 1, 0},
                        { 0 ,0 ,0, 0}
                                });



                Matrix<double> Tt = T.Transpose(); //transpose
                K_eG = ke.Multiply(T);  //fikset denne


                K_eG = Tt.Multiply(K_eG);
                K_eG = mat * K_eG;  //global element stivehetsmatrise
                ke = mat * ke;      //lokal element stivhetsmatrise



                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                //upper left corner of k-matrix
                K_tot[node1 * 2, node1 * 2] += K_eG[0, 0];
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

        private static void CreateReducedStiffnesMatrix(List<Point3d> _pts ,List<int> _BC, Matrix<double> K_tot, Vector<double> R, out Matrix<double> K_red, out Vector<double> R_red, out List<double> ptsK )
        {

            
            int sum = 0;
            List<int> index = new List<int>();


            int num = K_tot.RowCount;
            

            List<int> keep = new List<int>();
            for (int i = 0; i < num; i++)
            {
                if (!_BC.Contains(i))
                {
                    keep.Add(i);
                }
            }



            
            List<double> ptsCoord = new List<double>();
            for (int i = 0; i<_pts.Count;i++)
            {
                ptsCoord.Add(_pts[i].X);
                ptsCoord.Add(_pts[i].Z);
            }

            List<double> ptsKeep = new List<double>();
            List<double> R_redu = new List<double>();
            Matrix<double> K_redu = SparseMatrix.OfArray(new double[keep.Count, keep.Count]);

            for (int i = 0; i < keep.Count  ; i++)
            {

                ptsKeep.Add(ptsCoord[keep[i]]);
                R_redu.Add(R[keep[i]]);
                for (int j = 0; j < keep.Count; j++)
                {

                    K_redu[i, j] = K_tot[keep[i], keep[j]];
                    

                }


            }

            /*
            Matrix<double> K_ = SparseMatrix.OfArray(new double[_BC.Count, _BC.Count]);
            K_ = K_tot.RemoveRow(0, 1, 4, 5).RemoveColumn(0, 1, 4, 5);

            
            for (int i = 0; i < _BC.Count; i++)
            {
                //int sizeOfMatrix = K_tot.GetLength(0);
                int blockedIndex = _BC[i];



                K_ = K_tot.RemoveRow(blockedIndex).RemoveColumn(blockedIndex);
                
                
                
                R.ToColumnMatrix().RemoveRow(blockedIndex);
                
                */



            /*

            for (int r = 0; r < sizeOfMatrix; r++)
            {
                K_tot[blockedIndex, r] = 0;
                K_tot[r, blockedIndex] = 0;


            }


            K_tot[blockedIndex, blockedIndex]=1;





        }
        */
            
            Vector<double> R_reduu = Vector<double>.Build.DenseOfArray(R_redu.ToArray());
            K_red = K_redu;
            R_red = R_reduu;
            ptsK = ptsKeep;
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