using Grasshopper;
using Grasshopper.Kernel;
using System;
using System.Drawing;

namespace Master
{
    public class MasterInfo : GH_AssemblyInfo
    {
        public override string Name => "Master";

        //Return a 24x24 pixel bitmap to represent this GHA library.
        public override Bitmap Icon => null;

        //Return a short string describing the purpose of this GHA library.
        public override string Description => "";

        public override Guid Id => new Guid("F58C6254-0B7B-4EA2-86ED-21ACA3EDD9D8");

        //Return a string identifying you or your company.
        public override string AuthorName => "";

        //Return a string representing your preferred contact details.
        public override string AuthorContact => "";
    }
}