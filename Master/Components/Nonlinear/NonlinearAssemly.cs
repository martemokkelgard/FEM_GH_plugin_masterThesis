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
using MathNet.Numerics.Integration;



namespace Master.Components
{
    public class NonlinearAssembly2 : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Assembly2DTruss class.
        /// </summary>
        public NonlinearAssembly2()
          : base("NonlinearAssembly", "Nickname",
              "Description",
              "Panda", "3DBeam")
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
            pManager.AddNumberParameter("x-position", "Pos", "x-Position", GH_ParamAccess.item);
            pManager.AddGenericParameter("ElementStiffnessMatrixClass", "keClass", "matrixClass elementstiffnessMatrix and degree of shapefunctions", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Displacement [mm] in x", "def", "Displacement [mm]", GH_ParamAccess.list);
            pManager.AddPointParameter("Rotation [rad] in x", "R", "Forces in points", GH_ParamAccess.list);
            pManager.AddPointParameter("Forces [N]", "F", "Forces in points", GH_ParamAccess.list);
            pManager.AddPointParameter("Moment [Nmm]", "R", "Rotation in points", GH_ParamAccess.list);
            pManager.AddPointParameter("DisplacementOfNodes", "displ", "Deformed geometry", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strain in x", "E", "strain ", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stress [N/mm^2] in x", "S", "stress [N/mm^2] ", GH_ParamAccess.list);
            pManager.AddNumberParameter("Moment in x", "M", "Moment ", GH_ParamAccess.item);
            pManager.AddNumberParameter("umax", "Umax", "maximum displacement ", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {


            //input


            List<BeamClassNonLin> beams = new List<BeamClassNonLin>();
            List<LoadClass> lc = new List<LoadClass>();
            List<BcClass> bcc = new List<BcClass>();
            List<NodeClass> nodes = new List<NodeClass>();
            MatrixClass mc = new MatrixClass();
            double x = new double();


            DA.GetDataList(0, beams);
            DA.GetDataList(1, bcc);
            DA.GetDataList(2, lc);
            DA.GetData(3, ref x);
            DA.GetData(4, ref mc);


            Matrix<double> ke = mc.elementMatrix;
            double deg = mc.deg;
            Matrix<double> N = mc.N;
            Matrix<double> dN = mc.dN;

            List<double> param = new List<double>();
            for (int i = 0; i < deg + 1; i++)
            {
                param.Add(i / deg);
            }

            List<Point3d> pts = new List<Point3d>();  //making pointList from lines

            for (int i = 0; i < beams.Count; i++)
            {
                Curve L1 = beams[i].axis;
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


            List<int> BCList = CreateBCList(bcc, pts);

            Vector<double> LoadList = CreateLoadList(lc, pts);



            CreateGlobalStiffnesMatrix(beams, ke, pts, out Matrix<double> k_tot, out List<Matrix<double>> k_eg);

            var R = Vector<double>.Build.DenseOfArray(LoadList.ToArray());

            CreateReducedStiffnesMatrix(BCList, k_tot, R, out Matrix<double> K_red);

            //Matrix<double> invK = K_red.Inverse();

            //var def = invK.Multiply(R_red);
            var displNodes = new List<Point3d>();

            Vector<double> def = K_red.Cholesky().Solve(R);     //r


            


            CreateForces(beams, pts, def, k_eg, N, dN, out Vector<double> forces, out Vector<double> moment, out Vector<double> strain, out Vector<double> stress, out List<Vector<double>> u, out List<double> M, out double umax);


            for (int i = 0; i < pts.Count; i++)
            {

                displNodes.Add(new Point3d(pts[i].X + def[6 * i] / 1000, pts[i].Y + def[6 * i + 1] / 1000, pts[i].Z + def[6 * i + 2] / 1000));
            }



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


                if (Math.Abs(moment[i]) <= 0.00001)
                {
                    moment[i] = 0;
                }

            }



            foreach (Vector uu in u)
            {

                int tol = 6;
                disp_lst.Add(new Point3d(Math.Round(uu[0], tol), Math.Round(uu[1], tol), Math.Round(uu[2], tol)));
                rot_lst.Add(new Point3d(Math.Round(uu[3], tol), Math.Round(uu[4], tol), Math.Round(uu[5], tol)));
            }


            // endrer output til liste med punkter
            for (int i = 0; i < pts.Count; i++)
            {


                if (BCList.Contains(i * 6) | BCList.Contains(i * 6 + 1) | BCList.Contains(i * 6 + 2) | BCList.Contains(i * 6 + 3) | BCList.Contains(i * 6 + 4) | BCList.Contains(i * 6 + 5))
                {
                    force_lst.Add(new Point3d(forces[i * 3], forces[i * 3 + 1], forces[i * 3 + 2]));
                    mom_lst.Add(new Point3d(moment[i * 3], moment[i * 3 + 1], moment[i * 3 + 2]));

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
            DA.SetDataList(6, def);
            DA.SetData(7, M);
            DA.SetData(8, umax);


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
            /*
            List<int> globalIds = new List<int>();
            foreach (var b in _Bars)
            {
                globalIds.Add(b.startNode.Id);
                globalIds.Add(b.midNode.Id);
                globalIds.Add(b.endNode.Id);
            }
            globalIds.Distinct();
            
            int sizeOfR = globalIds.Count*6;
            */



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

        private static void CreateForces(List<BeamClassNonLin> bars, List<Point3d> points, Vector<double> _def, List<Matrix<double>> k_eg, Matrix<double> N, Matrix<double> dN, out Vector<double> forces, out Vector<double> moment, out Vector<double> strain, out Vector<double> stress, out List<Vector<double>> _u, out List<double> M, out double umax)
        {
            //Matrix<double> k_eG = DenseMatrix.OfArray(new double[6, 6]);
            Vector<double> u = SparseVector.OfEnumerable(new double[6]);
            Vector<double> eps = SparseVector.OfEnumerable(new double[6]);
            Vector<double> sigma = SparseVector.OfEnumerable(new double[6]);
            Vector<double> S;

            Vector<double> v = SparseVector.OfEnumerable(new double[18]);
            var _M = new double();

            Vector<double> mom = SparseVector.OfEnumerable(new double[points.Count * 3]);
            Vector<double> forc = SparseVector.OfEnumerable(new double[points.Count * 3]);

            double u_Max = 0;


            List<Vector<double>> u_lst = new List<Vector<double>>();
            List<double> M_lst = new List<double>();


            foreach (BeamClassNonLin b in bars)
            {


                var E = b.material.youngsModolus;
                var my = b.material.v;
                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;
                int node3 = b.midNodes.Id;

                v[0] = _def[node1 * 6];
                v[1] = _def[node1 * 6 + 1];
                v[2] = _def[node1 * 6 + 2];
                v[3] = _def[node1 * 6 + 3];
                v[4] = _def[node1 * 6 + 4];
                v[5] = _def[node1 * 6 + 5];

                v[6] = _def[node3 * 6];
                v[7] = _def[node3 * 6 + 1];
                v[8] = _def[node3 * 6 + 2];
                v[9] = _def[node3 * 6 + 3];
                v[10] = _def[node3 * 6 + 4];
                v[11] = _def[node3 * 6 + 5];

                v[12] = _def[node2 * 6];
                v[13] = _def[node2 * 6 + 1];
                v[14] = _def[node2 * 6 + 2];
                v[15] = _def[node2 * 6 + 3];
                v[16] = _def[node2 * 6 + 4];
                v[17] = _def[node2 * 6 + 5];


                //displacement

                var B = dN;
                u = N.Multiply(v);


                u_lst.Add(u);

                var u_max = new double();
                u_max = Math.Sqrt(Math.Pow(u[0], 2) + Math.Pow(u[1], 2) + Math.Pow(u[2], 2));

                if (u_max > u_Max)
                {
                    u_Max = u_max;
                }

                eps = B.Multiply(v);        //strain (3)

                sigma = E * eps;            //stress = E*3
                S = k_eg[b.Id].Multiply(v);

                _M = E * b.section.Iy * eps[2];

                M_lst.Add(_M);


                forc[node1 * 3] += S[0];
                forc[node1 * 3 + 1] += S[1];
                forc[node1 * 3 + 2] += S[2];

                forc[node3 * 3] += S[6];
                forc[node3 * 3 + 1] += S[7];
                forc[node3 * 3 + 2] += S[8];

                forc[node2 * 3] += S[12];
                forc[node2 * 3 + 1] += S[13];
                forc[node2 * 3 + 2] += S[14];


                mom[node1 * 3] += S[3];
                mom[node1 * 3 + 1] += S[4];
                mom[node1 * 3 + 2] += S[5];

                mom[node3 * 3] += S[9];
                mom[node3 * 3 + 1] += S[10];
                mom[node3 * 3 + 2] += S[11];


                mom[node2 * 3] += S[15];
                mom[node2 * 3 + 1] += S[16];
                mom[node2 * 3 + 2] += S[17];



            }


            forces = forc;
            moment = mom;
            strain = eps;
            stress = sigma;
            umax = u_Max;
            _u = u_lst;
            M = M_lst;

        }


        private static void CreateGenerelizedShapeFunc(List<BeamClassNonLin> bars, double X, out Matrix<double> N, out Matrix<double> dN)
        {

            Matrix<double> _N = DenseMatrix.OfArray(new double[6, 18]);
            Matrix<double> _dN = DenseMatrix.OfArray(new double[6, 18]);



            foreach (BeamClassNonLin b in bars)
            {


                Curve currentLine = b.axis;
                double L = currentLine.GetLength() * 1000.00;
                var x = X * L;


                double N1 = 2.00 * Math.Pow(x, 2) / (Math.Pow(L, 2)) - 3.00 * x / L + 1.00;                                //axial translation in node 1
                double N2 = 4.00 * x / L - 4.00 * Math.Pow(x, 2) / (Math.Pow(L, 2));                                          //axial translation in midside node
                double N3 = 2.00 * Math.Pow(x, 2) / (Math.Pow(L, 2)) - x / L;                                            //axial translation in node 2

                double N4 = 66.00 * Math.Pow(x, 3) / (Math.Pow(L, 3)) - 23.00 * Math.Pow(x, 2) / (Math.Pow(L, 2)) - 68.00 * Math.Pow(x, 4) / (Math.Pow(L, 4)) + 24.00 * Math.Pow(x, 5) / (Math.Pow(L, 5)) + 1.00;             //translation node 1
                double N5 = -x + 6.00 * Math.Pow(x, 2) / L - 13.00 * Math.Pow(x, 3) / (Math.Pow(L, 2)) + 12.00 * Math.Pow(x, 4) / (Math.Pow(L, 3)) - 4.00 * Math.Pow(x, 5) / (Math.Pow(L, 4));                                              //rotation node 1

                double N6 = 16.00 * Math.Pow(x, 2) / (Math.Pow(L, 2)) - 32.00 * Math.Pow(x, 3) / (Math.Pow(L, 3)) + 16.00 * Math.Pow(x, 4) / (Math.Pow(L, 4));                                                  //translation midside node
                double N7 = -32.00 * Math.Pow(x, 3) / Math.Pow(L, 2) + 8.00 * Math.Pow(x, 2) / L + 40.00 * Math.Pow(x, 4) / (Math.Pow(L, 3)) - 16.00 * Math.Pow(x, 5) / (Math.Pow(L, 4));                           //rotation midside node

                double N8 = 7.00 * Math.Pow(x, 2) / (Math.Pow(L, 2)) - 34.00 * Math.Pow(x, 3) / (Math.Pow(L, 3)) + 52.00 * Math.Pow(x, 4) / (Math.Pow(L, 4)) - 24.00 * Math.Pow(x, 5) / (Math.Pow(L, 5));          //translation node 2
                double N9 = -5.00 * Math.Pow(x, 3) / (Math.Pow(L, 2)) + Math.Pow(x, 2) / L + 8.00 * Math.Pow(x, 4) / (Math.Pow(L, 3)) - 4.00 * Math.Pow(x, 5) / (Math.Pow(L, 4));                                 //rotation node 2

                double dN1 = 4.00 * x / Math.Pow(L, 2) - 3.00 / L;
                double dN2 = 4.00 / L - 8.00 * x / Math.Pow(L, 2);
                double dN3 = 4.00 * x / Math.Pow(L, 2) - 1.00 / L;

                double dN4 = 198.00 * Math.Pow(x, 2) / Math.Pow(L, 3) - 46.00 * x / Math.Pow(L, 2) - 272.00 * Math.Pow(x, 3) / Math.Pow(L, 4) + 120.00 * Math.Pow(x, 4) / Math.Pow(L, 5);
                double dN5 = -39.00 * Math.Pow(x, 2) / Math.Pow(L, 2) + 12.00 * x / L + 48.00 * Math.Pow(x, 3) / Math.Pow(L, 3) - 20.00 * Math.Pow(x, 4) / Math.Pow(L, 4) - 1.00;
                double dN6 = 32.00 * x / Math.Pow(L, 2) - 96.00 * Math.Pow(x, 2) / Math.Pow(L, 3) + 64.00 * Math.Pow(x, 3) / Math.Pow(L, 4);
                double dN7 = -96.00 * Math.Pow(x, 2) / Math.Pow(L, 2) + 16.00 * x / L + 160.00 * Math.Pow(x, 3) / Math.Pow(L, 3) - 80.00 * Math.Pow(x, 4) / Math.Pow(L, 4);
                double dN8 = 14.00 * x / Math.Pow(L, 2) - 102.00 * Math.Pow(x, 2) / Math.Pow(L, 3) + 208.00 * Math.Pow(x, 3) / Math.Pow(L, 4) - 120.00 * Math.Pow(x, 4) / Math.Pow(L, 5);
                double dN9 = -15.00 * Math.Pow(x, 2) / Math.Pow(L, 2) + 2.00 * x / L + 32.00 * Math.Pow(x, 3) / Math.Pow(L, 3) - 20.00 * Math.Pow(x, 4) / Math.Pow(L, 4);

                double ddN4 = 396.00 * x / Math.Pow(L, 3) - 46.00 / Math.Pow(L, 2) - 816.00 * Math.Pow(x, 2) / Math.Pow(L, 4) + 480.00 * Math.Pow(x, 3) / Math.Pow(L, 5);
                double ddN5 = -78.00 * x / Math.Pow(L, 2) + 12.00 / L + 144.00 * Math.Pow(x, 2) / Math.Pow(L, 3) - 80.00 * Math.Pow(x, 3) / Math.Pow(L, 4);

                double ddN6 = 32.00 / Math.Pow(L, 2) - 192.00 * x / Math.Pow(L, 3) + 192.00 * Math.Pow(x, 2) / Math.Pow(L, 4);
                double ddN7 = -192.00 * x / Math.Pow(L, 2) + 16.00 / L + 480.00 * Math.Pow(x, 2) / Math.Pow(L, 3) - 320.00 * Math.Pow(x, 3) / Math.Pow(L, 4);

                double ddN8 = 14.00 / Math.Pow(L, 2) - 204.00 * x / Math.Pow(L, 3) + 624.00 * Math.Pow(x, 2) / Math.Pow(L, 4) - 480.00 * Math.Pow(x, 3) / Math.Pow(L, 5);
                double ddN9 = -30.00 * x / Math.Pow(L, 2) + 2.00 / L + 96.00 * Math.Pow(x, 2) / Math.Pow(L, 3) - 80.00 * Math.Pow(x, 3) / Math.Pow(L, 4);

                // første rad blir ux = N1*ux1 + N2*ux2
                //andre rad blir uy = N3*uy1 + N5*uy2 etc.

                _N = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                {N1,   0,      0,     0,     0,      0,    N2,      0,     0,     0,     0,     0,    N3,    0,     0,    0,     0,      0},
                {0,   N4,      0,     0,     0,     N5,     0,     N6,     0,     0,     0,    N7,    0,    N8,     0,    0,     0,     N9},
                {0,    0,     N4,     0,    -N5,     0,     0,     0,     N6,     0,   -N7,     0,    0,     0,     N8,   0,    -N9,     0},
                {0,    0,      0,    N1,     0,      0,     0,     0,      0,    N2,     0,     0,    N3,    0,     0,    0,     0,      0},
                {0,    dN4,    0,     0,     0,    dN5,     0,     dN6,    0,     0,     0,   dN7,    0,    dN8,    0,    0,     0,    dN9},
                {0,    0,    dN4,     0,   -dN5,     0,     0,     0,    dN6,     0,   -dN7,    0,    0,     0,    dN8,   0,    -dN9,    0},

                });



                //N1 and N2 is linear shape function and is derived one (Translationin x dir) while the rest is derived twice for beams. 

                _dN = Matrix<double>.Build.DenseOfArray(new double[,]
                {
                {dN1,   0,     0,     0,      0,      0,    dN2,     0,     0,      0,     0,      0,     dN3,     0,     0,     0,      0,       0},
                {0,   dN4,     0,     0,      0,    dN5,     0,    dN6,     0,      0,     0,     dN7,     0,     dN8,    0,     0,      0,     dN9},
                {0,    0,    dN4,     0,   -dN5,      0,     0,      0,    dN6,     0,   -dN7,     0,      0,      0,    dN8,    0,    -dN9,      0},
                {0,    0,     0,    dN1,      0,      0,     0,      0,     0,    dN2,     0,      0,      0,      0,     0,    dN3,     0,       0},
                {0,    ddN4,    0,    0,      0,    ddN5,    0,     ddN6,   0,      0,     0,    ddN7,     0,    ddN8,    0,     0,      0,     ddN9},
                {0,     0,    ddN4,   0,   -ddN5,     0,     0,      0,    ddN6,    0,   -ddN7,    0,      0,      0,    ddN8,   0,    -ddN9,     0},

                });



            }

            N = _N;
            dN = _dN;
        }






        private static void CreateGlobalStiffnesMatrix(List<BeamClassNonLin> bars, Matrix<double> _ke, List<Point3d> points, out Matrix<double> k_tot, out List<Matrix<double>> k_eg)

        {

            int dofs = (points.Count) * 6;
            Matrix<double> K_tot = DenseMatrix.OfArray(new double[dofs, dofs]);
            Matrix<double> K_eG = DenseMatrix.OfArray(new double[18, 18]);

            List<Matrix<double>> ke_lst = new List<Matrix<double>>();


            foreach (BeamClassNonLin b in bars)
            {


                Curve currentLine = b.axis;
                double L = currentLine.GetLength() * 1000.00;



                double X1 = 7.00 / 3.00 * b.section.CSA * b.material.youngsModolus / L;
                double X2 = -8.00 / 3.00 * b.section.CSA * b.material.youngsModolus / L;
                double X3 = 1.00 / 3.00 * b.section.CSA * b.material.youngsModolus / L;
                double X4 = 16.00 / 3.00 * b.section.CSA * b.material.youngsModolus / L;


                double Y1 = (5092.00 / 35.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 3));
                double Y2 = (-1138.00 / 35.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 2));
                double Y3 = (-512.00 / 5.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 3));
                double Y4 = (-384.00 / 7.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 2));
                double Y5 = (-1508.00 / 35.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 3));
                double Y6 = (-242.00 / 35.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 2));
                double Y7 = (332.00 / 35.00) * b.material.youngsModolus * b.section.Iz / L;
                double Y8 = (128.00 / 5.00) * b.material.youngsModolus * b.section.Iz / (Math.Pow(L, 2));
                double Y9 = (64.00 / 7.00) * b.material.youngsModolus * b.section.Iz / L;
                double Y10 = (38.00 / 35.00) * b.material.youngsModolus * b.section.Iz / L;
                double Y11 = (1024.00 / 5.00) * b.material.youngsModolus * b.section.Iz / Math.Pow(L, 3);
                double Y12 = (256.00 / 7.00) * b.material.youngsModolus * b.section.Iz / L;

                double Z1 = (5092.00 / 35.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 3));
                double Z2 = (-1138.00 / 35.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 2));
                double Z3 = (-512.00 / 5.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 3));
                double Z4 = (-384.00 / 7.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 2));
                double Z5 = (-1508.00 / 35.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 3));
                double Z6 = (-242.00 / 35.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 2));
                double Z7 = (332.00 / 35.00) * b.material.youngsModolus * b.section.Iy / L;
                double Z8 = (128.00 / 5.00) * b.material.youngsModolus * b.section.Iy / (Math.Pow(L, 2));
                double Z9 = (64.00 / 7.00) * b.material.youngsModolus * b.section.Iy / L;
                double Z10 = (38.00 / 35.00) * b.material.youngsModolus * b.section.Iy / L;
                double Z11 = (1024.00 / 5.00) * b.material.youngsModolus * b.section.Iy / Math.Pow(L, 3);
                double Z12 = (256.00 / 7.00) * b.material.youngsModolus * b.section.Iy / L;

                double C1 = (7.00 / 3.00) * (b.material.G * b.section.J) / L;
                double C2 = (-8.00 / 3.00) * (b.material.G * b.section.J) / L;
                double C3 = (1.00 / 3.00) * (b.material.G * b.section.J) / L;
                double C4 = (16.00 / 3.00) * (b.material.G * b.section.J) / L;

                Point3d p1 = b.startNode.pt;
                Point3d p2 = b.endNode.pt;
                Point3d p3 = b.midNode.pt;


                double xl = (p2.X - p1.X);
                double yl = (p2.Y - p1.Y);
                double zl = (p2.Z - p1.Z);

                double l = currentLine.GetLength();
                double den = l * Math.Pow(Math.Pow(xl, 2) + Math.Pow(yl, 2), 0.5);

                double cx = xl / l;
                double cy = yl / l;
                double cz = zl / l;

                double s = (p2.Z - p1.Z) / l;
                double c = (Math.Pow(Math.Pow((p2.X - p1.X), 2) + Math.Pow((p2.Y - p1.Y), 2), 0.5)) / l;

                Matrix<double> t = DenseMatrix.OfArray(new double[,]
                {
                        {cx,                                     cy,                   cz},
                        {-(xl*zl*s + l*yl*c) / den,     -(yl*zl*s - l*xl*c) / den,     den*s/(l*l)},
                        {(xl*zl*c - l*yl*s)/den,         (yl*zl*c + l*xl*s) / den,       -den*c / (l*l)},
                });

                var T_t = t.DiagonalStack(t);
                var T = T_t.DiagonalStack(T_t);
                T = T.DiagonalStack(T_t);



                Matrix<double> ke = DenseMatrix.OfArray(new double[,]
                {
                        {X1,  0,   0,   0,    0,    0,   X2,    0,    0,   0,    0,    0,    X3,    0,     0,    0,     0,     0},
                        {0,  Y1,   0,   0,    0,   Y2,    0,   Y3,    0,   0,    0,    Y4,    0,    Y5,    0,    0,     0,    Y6},
                        {0,   0,   Z1,  0,   Z2,    0,    0,    0,   Z3,   0,    Z4,   0,     0,    0,     Z5,   0,    Z6,    0},
                        {0,   0,   0,   C1,   0,    0,    0,    0,   0,   C2,    0,    0,     0,    0,     0,    C3,    0,    0},
                        {0,   0,  Z2,   0,   Z7,    0,    0,    0,   Z8,   0,    Z9,   0,     0,    0,    -Z6,   0,    Z10,   0},
                        {0,   Y2,  0,   0,    0,   Y7,    0,   Y8,    0,   0,    0,    Y9,    0,   -Y6,    0,    0,     0,   Y10},
                        {X2,  0,   0,   0,    0,    0,   X4,    0,    0,   0,    0,    0,    X2,    0,     0,    0,     0,     0},
                        {0,  Y3,   0,   0,    0,    Y8,   0,   Y11,   0,   0,    0,    0,     0,   Y3,     0,    0,     0,   -Y8},
                        {0,   0,  Z3,   0,   Z8,    0,    0,    0,   Z11,  0,    0,    0,     0,    0,    Z3,    0,    -Z8,   0},
                        {0,   0,   0,  C2,    0,    0,    0,    0,    0,   C4,   0,    0,     0,    0,     0,    C2,    0,     0},
                        {0,   0,  Z4,   0,   Z9,    0,    0,    0,    0,   0,   Z12,   0,     0,    0,    -Z4,   0,    Z9,    0},
                        {0,  Y4,   0,   0,    0,   Y9,    0,    0,    0,   0,    0,    Y12,   0,   -Y4,    0,    0,     0,    Y9},
                        {X3,  0,   0,   0,    0,    0,   X2,    0,    0,   0,    0,    0,    X1,    0,     0,    0,     0,     0},
                        {0,  Y5,   0,   0,    0,   -Y6,   0,    Y3,   0,   0,    0,   -Y4,    0,   Y1,     0,    0,     0,   -Y2},
                        {0,   0,  Z5,   0,   -Z6,   0,    0,    0,   Z3,   0,   -Z4,   0,     0,    0,     Z1,   0,    -Z2,   0},
                        {0,   0,   0,  C3,    0,    0,    0,    0,    0,   C2,   0,    0,     0,    0,     0,   C1,     0,     0},
                        {0,   0,  Z6,   0,   Z10,   0,    0,    0,   -Z8,  0,   Z9,    0,     0,    0,    -Z2,   0,    Z7,    0},
                        {0,  Y6,   0,   0,    0,   Y10,   0,  -Y8,    0,   0,    0,    Y9,    0,   -Y2,    0,    0,     0,    Y7},


                });
                Matrix<double> Tt = T.Transpose(); //transpose
                Matrix<double> KG = Tt.Multiply(ke);
                K_eG = KG.Multiply(T);

                ke_lst.Add(K_eG);

                int node1 = b.startNode.Id;
                int node2 = b.endNode.Id;
                int node3 = b.midNode.Id;

                int tol = 15;

                for (int i = 0; i < K_eG.RowCount / 3; i++)
                {
                    for (int j = 0; j < K_eG.ColumnCount / 3; j++)
                    {
                        K_tot[node1 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i, j], tol);
                        K_tot[node1 * 6 + i, node3 * 6 + j] += Math.Round(K_eG[i, j + 6], tol);
                        K_tot[node1 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i, j + 12], tol);

                        K_tot[node3 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i + 6, j], tol);
                        K_tot[node3 * 6 + i, node3 * 6 + j] += Math.Round(K_eG[i + 6, j + 6], tol);
                        K_tot[node3 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i + 6, j + 12], tol);

                        K_tot[node2 * 6 + i, node1 * 6 + j] += Math.Round(K_eG[i + 12, j], tol);
                        K_tot[node2 * 6 + i, node3 * 6 + j] += Math.Round(K_eG[i + 12, j + 6], tol);
                        K_tot[node2 * 6 + i, node2 * 6 + j] += Math.Round(K_eG[i + 12, j + 12], tol);

                    }

                }


            }

            k_tot = K_tot;
            k_eg = ke_lst;

        }

        private static void CreateReducedStiffnesMatrix(List<int> _BC, Matrix<double> K_tott, Vector<double> Rvec, out Matrix<double> K_red)
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
            get { return new Guid("A4B9D9D7-4125-4BAA-9AA1-4CE34A7471CF"); }
        }
    }
}