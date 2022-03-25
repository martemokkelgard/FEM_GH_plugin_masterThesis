﻿using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Linq;
using MathNet.Numerics.LinearAlgebra.Factorization;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using MathNet.Numerics.Integration;



namespace Master.Components
{
    public class NonlinearAssembly : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Assembly2DTruss class.
        /// </summary>
        public NonlinearAssembly()
          : base("NonlinearAssembly", "Nickname",
              "Description",
              "Løve", "3DBeam")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Bar", "B", "BarClass object", GH_ParamAccess.list);
            pManager.AddGenericParameter("Boundary Conditions", "BC", "BcClass object", GH_ParamAccess.list);
            pManager.AddGenericParameter("Loads", "L", "LoadClass object", GH_ParamAccess.list);
            pManager.AddMatrixParameter("Kematrix","Ke","Ke_matrix from MathCad",GH_ParamAccess.list);
            pManager.AddMatrixParameter("Nmatrix", "N", "N_matrix from MathCad", GH_ParamAccess.item);
            pManager.AddMatrixParameter("dNmatrix", "dN", "dN_matrix from MathCad", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Displacement [mm]", "def", "Displacement [mm]", GH_ParamAccess.list);
            pManager.AddPointParameter("Rotation [rad]", "R", "Forces in points", GH_ParamAccess.list);
            pManager.AddPointParameter("Forces [N]", "F", "Forces in points", GH_ParamAccess.list);
            pManager.AddPointParameter("Moment [Nmm]", "R", "Rotation in points", GH_ParamAccess.list);
            pManager.AddPointParameter("DisplacementOfNodes", "displ", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain", "E", "strain ", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stress [N/mm^2] ", "S", "stress [N/mm^2] ", GH_ParamAccess.list);
            pManager.AddNumberParameter("displ ", "d", " displ. from form function ", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {


            //input


            List<BarClass> bars = new List<BarClass>();
            List<LoadClass> lc = new List<LoadClass>();
            List<BcClass> bcc = new List<BcClass>();
            //created in a differen component from MathCad
            Matrix<double> ke = DenseMatrix.OfArray(new double[12, 12]);
            Matrix<double> N = DenseMatrix.OfArray(new double[6, 12]);
            Matrix<double> dN = DenseMatrix.OfArray(new double[6, 12]);



            DA.GetDataList(0, bars);
            DA.GetDataList(1, bcc);
            DA.GetDataList(2, lc);
            DA.GetData(3, ref ke);      
            DA.GetData(4, ref N);  
            DA.GetData(5, ref dN);


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

            Vector<double> LoadList = CreateLoadList(lc, pts);


            CreateGlobalStiffnesMatrix(bars, pts, ke, out Matrix<double> k_tot, out Matrix<double> k_eg);

            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(BCList, k_tot, R, out Matrix<double> K_red, out Vector<double> R_red);

            Matrix<double> invK = K_red.Inverse();

            //var def = invK.Multiply(R_red);
            var displNodes = new List<Point3d>();

            Vector<double> def = K_red.Cholesky().Solve(R);     //r


            CreateForces(bars, pts, def, k_eg, N, dN, out Vector<double> forces, out Vector<double> rotation, out Vector<double> strain, out Vector<double> stress, out List<Vector<double>> displ);


            var displace = DataTreeFromVectorList(displ);

            for (int i = 0; i < pts.Count; i++)
            {

                displNodes.Add(new Point3d(pts[i].X + def[3 * i] / 1000, pts[i].Y + def[3 * i + 1] / 1000, pts[i].Z + def[3 * i + 2] / 1000));
            }

            //var outTree = DataTreeFromVectorList(forces);

            /*
            List<Line> liness = new List<Line>();

            for (int i = 0; i < displNodes.Count - 1; i++)
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
            */

            //lage lister med displacement

            List<Point3d> force_lst = new List<Point3d>();
            List<Point3d> mom_lst = new List<Point3d>();
            List<Point3d> disp_lst = new List<Point3d>();
            List<Point3d> rot_lst = new List<Point3d>();

            //gjør om små tall til 0
            for (int i = 0; i < forces.Count; i++)
            {
                if (Math.Abs(forces[i]) <= 0.00001)
                {
                    forces[i] = 0;
                }


                if (Math.Abs(rotation[i]) <= 0.00001)
                {
                    rotation[i] = 0;
                }

            }

            // endrer output til liste med punkter
            for (int i = 0; i < pts.Count; i++)
            {
                disp_lst.Add(new Point3d(def[i * 6], def[i * 6 + 1], def[i * 6 + 2]));
                rot_lst.Add(new Point3d(def[i * 6 + 3], def[i * 6 + 4], def[i * 6 + 5]));

                if (BCList.Contains(i * 6) | BCList.Contains(i * 6 + 1) | BCList.Contains(i * 6 + 2) | BCList.Contains(i * 6 + 3) | BCList.Contains(i * 6 + 4) | BCList.Contains(i * 6 + 5))
                {
                    force_lst.Add(new Point3d(forces[i * 3], forces[i * 3 + 1], forces[i * 3 + 2]));
                    mom_lst.Add(new Point3d(rotation[i * 3], rotation[i * 3 + 1], rotation[i * 3 + 2]));

                }


            }



            //output
            DA.SetDataList(0, disp_lst);
            
            DA.SetDataList(1, rot_lst);
            DA.SetDataList(2, force_lst);
            DA.SetDataList(3, mom_lst);
            //DA.SetDataTree(1, outTree);
            DA.SetDataList(4, displNodes);
            DA.SetDataList(5, strain);
            DA.SetDataList(6, stress);
            DA.SetDataTree(7, displace);


        }

        private List<int> CreateBCList(List<BcClass> _BcValue, List<Point3d> _Pts)  //making a list with indexes of fixed BC
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
                            BCsIndex.Add(6 * i);

                        if (b.uy)
                            BCsIndex.Add(6 * i + 1);

                        if (b.uz)
                            BCsIndex.Add(6 * i + 2);

                        if (b.rx)
                            BCsIndex.Add(6 * i + 3);

                        if (b.ry)
                            BCsIndex.Add(6 * i + 4);

                        if (b.rz)
                            BCsIndex.Add(6 * i + 5);
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

        private Vector<double> CreateLoadList(List<LoadClass> _lc, List<Point3d> _Pts)
        {

            Vector<double> LoadValue = SparseVector.OfEnumerable(new double[_Pts.Count * 6]);

            foreach (var load in _lc)

            {
                //List<Vector3d> LoadVec = new List<Vector3d>();

                //Vector<double> LoadValue = new Vector<double>();
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
                            if (load.Id1)
                            {
                                LoadValue[i * 6] += LoadVec.X;
                                LoadValue[i * 6 + 1] += LoadVec.Y;
                                LoadValue[i * 6 + 2] += LoadVec.Z;
                            }

                            else
                            {
                                LoadValue[i * 6 + 3] += load.mVal[0];
                                LoadValue[i * 6 + 4] += load.mVal[1];
                                LoadValue[i * 6 + 5] += load.mVal[2];
                            }

                        }
                        else
                        {
                            LoadValue[i * 6] += 0;
                            LoadValue[i * 6 + 1] += 0;
                            LoadValue[i * 6 + 2] += 0;
                            LoadValue[i * 6 + 3] += 0;
                            LoadValue[i * 6 + 4] += 0;
                            LoadValue[i * 6 + 5] += 0;
                        }

                    }

                }

            }

            return LoadValue;

        }

        private static void CreateForces(List<BarClass> bars, List<Point3d> points, Vector<double> _def, Matrix<double> k_eg, Matrix<double> N, Matrix<double> dN, out Vector<double> forces, out Vector<double> rotation, out Vector<double> strain, out Vector<double> stress, out List<Vector<double>> displ)
        {
            //Matrix<double> k_eG = DenseMatrix.OfArray(new double[6, 6]);
            Vector<double> u = SparseVector.OfEnumerable(new double[6]);
            Vector<double> eps = SparseVector.OfEnumerable(new double[3]);
            Vector<double> sigma = SparseVector.OfEnumerable(new double[3]);
            Vector<double> S;

            Vector<double> v = SparseVector.OfEnumerable(new double[12]);
            
            List<Vector<double>> u_lst = new List<Vector<double>>();

            Vector<double> rot = SparseVector.OfEnumerable(new double[points.Count * 3]);
            Vector<double> disp = SparseVector.OfEnumerable(new double[points.Count * 3]);
            

            foreach (BarClass b in bars)
            {
                
                var B = dN;
                var E = b.material.youngsModolus;
                var my = b.material.v;
                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                v[0] = _def[node1 * 6];
                v[1] = _def[node1 * 6 + 1];
                v[2] = _def[node1 * 6 + 2];
                v[3] = _def[node1 * 6 + 3];
                v[4] = _def[node1 * 6 + 4];
                v[5] = _def[node1 * 6 + 5];
                v[6] = _def[node2 * 6];
                v[7] = _def[node2 * 6 + 1];
                v[8] = _def[node2 * 6 + 2];
                v[9] = _def[node2 * 6 + 3];
                v[10] = _def[node2 * 6 + 4];
                v[11] = _def[node2 * 6 + 5];


                //displacement
                u = N.Multiply(v);

                u_lst.Add(u);

                eps = B.Multiply(v);

                var C = Matrix<double>.Build.DenseOfArray(new double[,]
                {

                    {1-my,   my,     my,      0,       0,      0},
                    {my,   1-my,     my,      0,       0,      0},
                    {my,    my,   (1-my),     0,       0,      0},
                    {0,      0,       0,  (1-2*my)/2,  0,      0},
                    {0,      0,       0,      0,  (1-2*my)/2,  0},
                    {0,      0,       0,      0,       0, (1-2*my)/2 },

                });

                C = E / ((1 + my) * (1 - 2 * my)) * C;

                sigma = C.Multiply(eps);

                S = k_eg.Multiply(v);

                disp[node1 * 3] += S[0];
                disp[node1 * 3 + 1] += S[1];
                disp[node1 * 3 + 2] += S[2];

                disp[node2 * 3] -= S[6];
                disp[node2 * 3 + 1] -= S[7];
                disp[node2 * 3 + 2] -= S[8];


                rot[node1 * 3] += S[3];
                rot[node1 * 3 + 1] += S[4];
                rot[node1 * 3 + 2] += S[5];

                rot[node2 * 3] -= S[9];
                rot[node2 * 3 + 1] -= S[10];
                rot[node2 * 3 + 2] -= S[11];




            }


            forces = disp;
            rotation = rot;
            strain = eps;
            stress = sigma;
            displ = u_lst;

        }




        private static void CreateGlobalStiffnesMatrix(List<BarClass> bars, List<Point3d> points, Matrix<double> ke, out Matrix<double> k_tot, out Matrix<double> k_eg)

        {

            int dofs = points.Count * 6;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[dofs, dofs]);
            Matrix<double> K_eG = DenseMatrix.OfArray(new double[6, 6]);
            

            foreach (BarClass b in bars)
            {


                Line currentLine = b.axis;
                double L = currentLine.Length * 1000;
                  


                Point3d p1 = new Point3d(Math.Round(currentLine.From.X, 6), Math.Round(currentLine.From.Y, 6), Math.Round(currentLine.From.Z, 6));
                Point3d p2 = new Point3d(Math.Round(currentLine.To.X, 6), Math.Round(currentLine.To.Y, 6), Math.Round(currentLine.To.Z, 6));


                double xl = (p2.X - p1.X);
                double yl = (p2.Y - p1.Y);
                double zl = (p2.Z - p1.Z);

                double l = currentLine.Length;
                double den = l * Math.Pow(Math.Pow(xl, 2) + Math.Pow(yl, 2), 0.5);

                double cx = xl / l;
                double cy = yl / l;
                double cz = zl / l;

                double s = (p2.Z - p1.Z) / l;
                double c = (Math.Pow(Math.Pow((p2.X - p1.X), 2) + Math.Pow((p2.Y - p2.Y), 2), 0.5)) / l;

                Matrix<double> t = DenseMatrix.OfArray(new double[,]
                {
                        {cx,                                     cy,                   cz},
                        {-(xl*zl*s + l*yl*c) / den,     -(yl*zl*s - l*xl*c) / den,     den*s/(l*l)},
                        {(xl*zl*c - l*yl*s)/den,         (yl*zl*c + l*xl*s) / den,       -den*c / (l*l)},
                });

                var T = t.DiagonalStack(t);
                T = T.DiagonalStack(T);

                
                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke);
                K_eG = KG.Multiply(T);


                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;

                for (int i = 0; i < K_eG.RowCount / 2; i++)
                {
                    for (int j = 0; j < K_eG.ColumnCount / 2; j++)
                    {
                        K_tot[node1 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i, j], 5);
                        K_tot[node1 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i, j + 6], 5);
                        K_tot[node2 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i + 6, j], 5);
                        K_tot[node2 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i + 6, j + 6], 5);

                    }

                }


            }

            k_tot = K_tot;
            k_eg = K_eG;

        }

        private static void CreateReducedStiffnesMatrix(List<int> _BC, Matrix<double> K_tott, Vector<double> Rvec, out Matrix<double> K_red, out Vector<double> R_red)
        {


            for (int i = 0; i < _BC.Count; i++)
            {
                int sizeOfMatrix = K_tott.ColumnCount;
                int blockedIndex = _BC[i];

                Rvec[blockedIndex] = 0;

                for (int r = 0; r < sizeOfMatrix; r++)
                {
                    K_tott[blockedIndex, r] = 0;
                    K_tott[r, blockedIndex] = 0;

                    if (K_tott[r, r] == 0)
                    {
                        K_tott[r, r] = 1;
                    }
                }
            }

            K_red = K_tott;
            R_red = Rvec;


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