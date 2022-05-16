using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

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
              "Panda", "2DTruss")
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

            pManager.AddPointParameter("Displacement [mm]", "def", "Displacement [mm] in nodes", GH_ParamAccess.list);
            pManager.AddPointParameter("Reaction Forces [N]", "F", "Forces in support points", GH_ParamAccess.list);
            pManager.AddPointParameter("Position of Supports", "PS", "Position for support points", GH_ParamAccess.list);
            pManager.AddPointParameter("Displacement of Nodes", "displ", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain", "S", "strain ", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stress [N/mm^2] ", "S", "stress [N/mm^2] ", GH_ParamAccess.list);


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
            List<BarClass2DTruss> bars = new List<BarClass2DTruss>();
            List<LoadClass2DTruss> lc = new List<LoadClass2DTruss>();
            List<BcClass2DTruss> bcc = new List<BcClass2DTruss>();



            DA.GetDataList(0, bars);
            DA.GetDataList(1, bcc);
            DA.GetDataList(2, lc);


            //code


            foreach (var load in lc)
            {
                LoadPts.Add(load.coordinate);
                LoadVec.Add(load.LoadVec);
            }


            List<Point3d> pts = new List<Point3d>();
            foreach (BarClass2DTruss b in bars)
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

            Vector<double> LoadList = CreateLoadList(lc, pts);
            List<int> BC = new List<int>();


            CreateGlobalStiffnesMatrix(bars, pts, out Matrix<double> k_tot);


            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(BCList, k_tot, out Matrix<double> K_red);




            var invK = K_red.Inverse();

            var def = invK.Multiply(R);     //r = K^-1 * R --> finding displacements

            var displNodes = new List<Point3d>();   //showing the displacement of the nodes in Rhino



            CreateForces(bars, pts, def, out Vector<double> forces);




            for (int i = 0; i < pts.Count; i++)
            {

                displNodes.Add(new Point3d(pts[i].X + def[2 * i]/1000 , pts[i].Y, pts[i].Z + def[2 * i + 1]/1000 ));
            }

            List<Line> liness = new List<Line>();

            
            for (int i = 0; i < displNodes.Count - 1; i++)
            {
                Line line = new Line(displNodes[i], displNodes[i + 1]);
                liness.Add(line);
            }
            
            List<double> strain = new List<double>();   //strain and stress
            List<double> stress = new List<double>();
            foreach (BarClass2DTruss b in bars)
            {
                double originLength = b.axis.Length;
                double deformedLength = liness[b.Id].Length;

                double dL = originLength - deformedLength;

                double e = dL / originLength;
                strain.Add(e);

                double s = e * b.material.youngsModolus;
                stress.Add(s);
            }

            //lage lister med displacement

            List<Point3d> force_lst = new List<Point3d>();

            List<Point3d> disp_lst = new List<Point3d>();


            List<Point3d> force_mom_pos = new List<Point3d>();

            //gjør om små tall til 0
            for (int i = 0; i < forces.Count; i++)
            {
                if (Math.Abs(forces[i]) <= 0.00001)
                {
                    forces[i] = 0;
                }
            }



            // endrer output til liste med punkter
            for (int i = 0; i < pts.Count; i++)
            {
                disp_lst.Add(new Point3d(def[i * 2], 0.00, def[i * 2 + 1]));


                if (BCList.Contains(i * 2) | BCList.Contains(i * 2 + 1))
                {
                    force_lst.Add(new Point3d(forces[i * 2], 0.00, forces[i * 2 + 1]));

                    force_mom_pos.Add(pts[i]);
                }

            }




            //output
            DA.SetDataList(0, disp_lst);

            DA.SetDataList(1, force_lst);

            DA.SetDataList(2, force_mom_pos);

            DA.SetDataList(3, displNodes);
            DA.SetDataList(4, strain);
            DA.SetDataList(5, stress);
        }

        /*
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
        */

        private List<int> CreateBCList(List<BcClass2DTruss> _BcValue, List<Point3d> _Pts)  //making a list with indexes of fixed BC
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
                        if (b.ux)
                            BCsIndex.Add(2 * i);

                        if (b.uz)
                            BCsIndex.Add(2 * i + 1);
                    }

                }

            }

            return BCsIndex;
        }

        private Vector<double> CreateLoadList(List<LoadClass2DTruss> _lc, List<Point3d> _Pts)
        {

            Vector<double> LoadValue = SparseVector.OfEnumerable(new double[_Pts.Count * 2]);

            foreach (var load in _lc)

            {

                List<Point3d> LoadPts = new List<Point3d>();

                if (load.Id == true)
                {
                    LoadPts.Add(load.coordinate);
                }


                else
                {
                    LoadPts.Add(load.stPt);
                    LoadPts.Add(load.enPt);
                }


                Vector3d LoadVec = load.LoadVec;

                for (int i = 0; i < _Pts.Count; i++)
                {
                    //this line was not so good
                    Point3d node = _Pts[i];


                    for (int j = 0; j < LoadPts.Count; j++)
                    {
                        Point3d loadLocation = LoadPts[j];
                        if (loadLocation.DistanceTo(node) < 0.00001)
                        {
                            LoadValue[i * 2] += LoadVec.X;
                            LoadValue[i * 2 + 1] += LoadVec.Z;

                        }
                        else
                        {
                            LoadValue[i * 2] += 0;
                            LoadValue[i * 2 + 1] += 0;

                        }

                    }

                }


            }

            return LoadValue;
        }

        private static void CreateForces(List<BarClass2DTruss> bars, List<Point3d> points, Vector<double> _def, out Vector<double> forces)
        {



            Vector<double> S;
            Vector<double> disp = SparseVector.OfEnumerable(new double[points.Count * 2]);



            foreach (BarClass2DTruss b in bars)
            {
                Vector<double> v = DenseVector.OfEnumerable(new double[4]);

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
                        { 0, 0, c, s},
                        { 0 ,0 ,-s, c}
                                    });

                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                                    {
                        { 1, 0, -1 ,0},
                        { 0, 0 ,0, 0},
                        { -1, 0, 1, 0},
                        { 0 ,0 ,0, 0}
                                    });


                ke = ke * mat;
                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke);
                Matrix<double> K_eG = KG.Multiply(T);

                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                v[0] = _def[node1 * 2];
                v[1] = _def[node1 * 2 + 1];
                v[2] = _def[node2 * 2];
                v[3] = _def[node2 * 2 + 1];

                S = K_eG.Multiply(v);

                disp[node1 * 2] += S[0];
                disp[node1 * 2 + 1] += S[1];

                disp[node2 * 2] -= S[2];
                disp[node2 * 2 + 1] -= S[3];


            }


            forces = disp;



        }

        private static void CreateGlobalStiffnesMatrix(List<BarClass2DTruss> bars, List<Point3d> points, out Matrix<double> k_tot)

        {

            int dofs = points.Count * 2;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[dofs, dofs]);
            //Matrix<double> K_eG = DenseMatrix.OfArray(new double[4, 4]);

            foreach (BarClass2DTruss b in bars)
            {


                Line currentLine = b.axis;
                double mat = (b.section.CSA * b.material.youngsModolus) / (currentLine.Length * 1000);

                Point3d p1 = currentLine.From;
                Point3d p2 = currentLine.To;

                double angle = Math.Atan2(p2.Z - p1.Z, p2.X - p1.X); 
                double c = Math.Cos(angle); 
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
                        { c,  s, 0, 0},
                        { -s, c, 0, 0},
                        { 0,  0, c, s},
                        { 0,  0,-s, c}
                                });

                Matrix<double> ke = SparseMatrix.OfArray(new double[,]
                                {
                        { 1, 0, -1 ,0},
                        { 0, 0 ,0, 0},
                        { -1, 0, 1, 0},
                        { 0 ,0 ,0, 0}
                                });



                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke * mat);
                Matrix<double> K_eG = KG.Multiply(T);


                //K_eG = mat * K_eG;  //global element stivehetsmatrise
                //ke = mat * ke;      //lokal element stivhetsmatrise

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


                K_tot[blockedIndex, blockedIndex] = 1;


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
            get { return new Guid("56CE183E-9717-44C8-8A80-88AA1991C6EE"); }
        }
    }
}