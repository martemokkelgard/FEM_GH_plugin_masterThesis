using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



namespace Master.Components
{
    public class Assembly2DBeam : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Assembly2DTruss class.
        /// </summary>
        public Assembly2DBeam()
          : base("Assembly2DBeam", "Nickname",
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
            pManager.AddNumberParameter("Defomation [mm]", "def", "Deformations of the bar in [mm]", GH_ParamAccess.list);
            pManager.AddNumberParameter("ReactionForces", "R", "Reaction Forces", GH_ParamAccess.tree);
            pManager.AddPointParameter("DisplacementOfNodes", "displ", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stress","S","stress of beam",GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain", "S", "strain oin beam", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            
            
            //input
            
            
            List<Point3d> LoadPts = new List<Point3d>();
            List<Vector3d> LoadVec = new List<Vector3d>();

            List<BarClass> bars = new List<BarClass>();
            List<LoadClass> lc = new List<LoadClass>();
            List<BcClass> bcc = new List<BcClass>();

            
            //int node1 = bars[1].endNode.Id;

            DA.GetDataList(0, bars);
            DA.GetDataList(1,  bcc);
            DA.GetDataList(2, lc);

          
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
                      

            Matrix<double> K_red;

            //Matrix<double> k_eG;
            Matrix<double> k_tot;
            //List<Vector<double>> forces = new List<Vector<double>>();
            



            CreateGlobalStiffnesMatrix(bars, pts, out k_tot);

            

            
            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(BCList, k_tot, out K_red);

            //Matrix<double> k_red = K_red.PointwiseRound();

            //Matrix<double> invK = DenseMatrix.OfArray(new double[9, 9]);
            Matrix<double> invK = K_red.Inverse();
            
            var def = invK.Multiply(R);

            var displNodes = new List<Point3d>();

            //Vector<double> def = K_red.Cholesky().Solve(R);


            CreateForces(bars, pts, def, out List<Vector<double>> forces);


            for (int i =0; i< pts.Count; i++)
            {

                
                displNodes.Add(new Point3d(pts[i].X + def[3 * i]/1000, pts[i].Y + def[3 * i + 1] / 1000, pts[i].Z + def[3 * i + 2]/1000));
            }

            var outTree = DataTreeFromVectorList(forces);

            
            List<Line> liness = new List<Line>();
            
            for (int i =0; i< displNodes.Count - 1; i++)
            {
                Line line = new Line(displNodes[i], displNodes[i + 1]);
                liness.Add(line);
            }
           
            List<double> strain = new List<double>();   //strain and stress
            List<double> stress = new List<double>();
            foreach (BarClass b in bars)
            {
                double originLength = b.axis.Length;
                double deformedLength = liness[b.Id].Length;

                double dL = originLength - deformedLength;

                double e = dL / originLength;
                strain.Add(e);

                double s = e * b.material.youngsModolus;
                stress.Add(s);
            }




            //output
            DA.SetDataList(0, def);
            DA.SetDataTree(1, outTree);
            DA.SetDataList(2, displNodes);
            DA.SetDataList(3, stress);
            DA.SetDataList(4, strain);
            

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
                            BCsIndex.Add(3*i);

                        if (b.uz)
                            BCsIndex.Add(3 * i + 1);

                        if (b.rx)
                            BCsIndex.Add(3*i+2);
                    }
                }

            }

            return BCsIndex;
        }

        private GH_Structure<GH_Number> DataTreeFromVectorList(List<Vector<double>> vecLst)
        {
            GH_Structure<GH_Number> tree = new GH_Structure<GH_Number>();

            int count = 0;
            foreach (Vector<double> vec in vecLst)
            {
                GH_Path path = new GH_Path(count);
                foreach (var num in vec.AsArray())
                {
                    tree.Append(new GH_Number(num), path);
                }

                count++;
            }

            return tree;
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
                        LoadValue.Add(_LoadValue[j][1]);
                        LoadValue.Add(_LoadValue[j][2]);
                    }
                    else
                    {
                        LoadValue.Add(0);
                        LoadValue.Add(0);
                        LoadValue.Add(0);
                    }

                }

            }

            return LoadValue;
        }

        private static void CreateForces(List<BarClass> bars, List<Point3d> points, Vector<double> _def, out List<Vector<double>> forces)
        {
            //Matrix<double> k_eG = DenseMatrix.OfArray(new double[6, 6]);
            Vector<double> S;
            List<Vector<double>> ST = new List<Vector<double>>();
            Vector<double> v = SparseVector.OfEnumerable(new double[6]);

            foreach (BarClass b in bars)
            {
                Line currentLine = b.axis;
                double L = currentLine.Length*1000;
                double LL = Math.Pow(L, 2);
                double mat = ( b.material.youngsModolus * b.section.I) / ( Math.Pow(L,3) );
                double my = (LL * b.section.CSA) / b.section.I;
                Point3d p1 = new Point3d(Math.Round(currentLine.From.X,5), Math.Round(currentLine.From.Y, 5), Math.Round(currentLine.From.Z, 5));
                Point3d p2 = new Point3d(Math.Round(currentLine.To.X, 5), Math.Round(currentLine.To.Y, 5), Math.Round(currentLine.To.Z, 5));
                //double lineLength = Math.Round(currentLine.Length, 6);


                // finding the cos value of the angle that projects the line to x,y,z axis (in 2D we use cos and sin of the same angle for x and z)

                double s = (p2.X - p1.X)/ currentLine.Length;
                double c = (p2.Z - p1.Z)/ currentLine.Length;
                

                Matrix<double> T = DenseMatrix.OfArray(new double[,]
                                    {
                        { c, s, 0, 0, 0 ,0},
                        {-s, c, 0, 0, 0 ,0},
                        {0, 0, 1, 0, 0 ,0 },
                        { 0, 0 ,0,c, s, 0},
                        { 0, 0 ,0,-s, c, 0},
                        { 0, 0 ,0, 0, 0, 1}
                                    });

                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                                    {
                        { my,0, 0, -my ,0, 0},
                        { 0, 12 , -6*L, 0,- 12, -6*L},
                        { 0, -6*L , 4*LL, 0, 6*L, 2*LL},
                        { -my, 0,0 ,my,0, 0},
                        { 0, -12 ,6*L, 0,12,6*L},
                        { 0 ,-6L ,2*LL, 0,6*L,4*LL}
                                    });
                ke = ke * mat;
                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke);
                Matrix<double> K_eG = KG.Multiply(T);


                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                v[0] = _def[node1 * 3];
                v[1] = _def[node1 * 3 +1];
                v[2] = _def[node1 * 3 + 2];
                v[3] = _def[node2 * 3];
                v[4] = _def[node2 * 3 +1];
                v[5] = _def[node2 * 3 + 2];

                S = K_eG.Multiply(v);
                ST.Add(S);

            }

            forces = new List<Vector<double>>(ST);
        }

        private static void CreateGlobalStiffnesMatrix(List<BarClass> bars, List<Point3d> points, out Matrix<double> k_tot)

        {

            int dofs = points.Count * 3;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[dofs, dofs]);
            //Matrix<double> K_eG = DenseMatrix.OfArray(new double[4, 4]);

            foreach (BarClass b in bars)
            {


                Line currentLine = b.axis;
                double L = currentLine.Length * 1000;
                double LL = Math.Pow(L, 2);
                double mat = (b.material.youngsModolus * b.section.I) / (Math.Pow(L, 3));
                double my = (LL * b.section.CSA) / b.section.I;
                Point3d p1 = new Point3d(Math.Round(currentLine.From.X, 5), Math.Round(currentLine.From.Y, 5), Math.Round(currentLine.From.Z, 5));
                Point3d p2 = new Point3d(Math.Round(currentLine.To.X, 5), Math.Round(currentLine.To.Y, 5), Math.Round(currentLine.To.Z, 5));
                double lineLength = Math.Round(currentLine.Length, 6);


                // finding the cos value of the angle that projects the line to x,y,z axis (in 2D we use cos and sin of the same angle for x and z)

                double s = (p2.X - p1.X) / lineLength;
                double c = (p2.Z - p1.Z) / lineLength;


                Matrix<double> T = SparseMatrix.OfArray(new double[,]
                                    {
                        { c, s, 0, 0, 0 ,0},
                        {-s, c, 0, 0, 0 ,0},
                        {0, 0, 1, 0, 0 ,0 },
                        { 0, 0 ,0,c, s, 0},
                        { 0, 0 ,0,-s, c, 0},
                        { 0, 0 ,0, 0, 0, 1}
                                    });

                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                                    {
                        { my,0, 0, -my ,0, 0},
                        { 0, 12 , -6*L, 0,- 12, -6*L},
                        { 0, -6*L , 4*LL, 0, 6*L, 2*LL},
                        { -my, 0,0 ,my,0, 0},
                        { 0, -12 ,6*L, 0,12,6*L},
                        { 0 ,-6L ,2*LL, 0,6*L,4*LL}
                                    });
                ke = ke * mat;
                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke);
                Matrix<double> K_eG = KG.Multiply(T);


                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                for (int i = 0; i < K_eG.RowCount / 2; i++)
                {
                    for (int j = 0; j < K_eG.ColumnCount / 2; j++)
                    {
                        K_tot[node1 * 3 + i, node1 * 3 + j] += K_eG[i, j];
                        K_tot[node1 * 3 + i, node2 * 3 + j] += K_eG[i, j+3];
                        K_tot[node2 * 3 + i, node1 * 3 + j] += K_eG[i+3, j];
                        K_tot[node2 * 3 + i, node2 * 3 + j] += K_eG[i+3, j+3];

                    }

                }


            }

            
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