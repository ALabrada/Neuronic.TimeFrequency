using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.TimeFrequency.Kernels;
using Neuronic.TimeFrequency.Testing.Properties;
using Neuronic.TimeFrequency.Wavelets;

namespace Neuronic.TimeFrequency.Testing
{
    [TestClass]
    public class KernelTests
    {
        [TestMethod]
        [DataRow(30d, "gen_kernel_cw_30", DisplayName = "Choi Williams sigma=30")]
        [DataRow(100d, "gen_kernel_cw_100", DisplayName = "Choi Williams sigma=100")]
        public void TestChoiWilliamsDistribution(double sigma, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            var expectedValues = new List<double[]>();
            using (var reader = new StringReader(valueList))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(c => (double)c).ToArray());
            }

            var kernel = new ChoiWilliamsDistribution {Sigma = sigma};
            var actualValues = kernel.Evaluate(expectedValues.Count);

            Assert.AreEqual(expectedValues.Count, actualValues.GetLength(0));
            for (int i = 0; i < expectedValues.Count; i++)
            {
                var row = new double[actualValues.GetLength(1)];
                Buffer.BlockCopy(actualValues, i * row.Length * sizeof(double), row, 0, row.Length * sizeof(double));
                Tools.AssertAreEqual(expectedValues[i], row, 1e-5);
            }
        }
    }
}
