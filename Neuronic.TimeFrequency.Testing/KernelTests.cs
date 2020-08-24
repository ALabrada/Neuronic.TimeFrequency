using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using Accord.IO;
using MathNet.Numerics;
using MathNet.Numerics.LinearAlgebra;
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
            var mat = new MatReader((byte[]) resources.GetObject(valueList));
            var expectedValues = CreateMatrix.DenseOfArray(mat.Read<double[,]>("g"));

            var kernel = new ChoiWilliamsDistribution {Sigma = sigma};
            var actualValues = CreateMatrix.DenseOfArray(kernel.Evaluate(expectedValues.RowCount));

            Assert.AreEqual(expectedValues.RowCount, actualValues.RowCount);
            Assert.AreEqual(expectedValues.ColumnCount, actualValues.ColumnCount);

            Assert.AreEqual(0, (expectedValues - actualValues).RowSums().Sum(), 1e-5);
        }

        [TestMethod]
        [DataRow("bart", 51, "gen_kernel_swvd_bart51_time", DisplayName = "Smoothed Wigner-Ville Distribution with bartlett window")]
        [DataRow("hann", 51, "gen_kernel_swvd_hann51_time", DisplayName = "Smoothed Wigner-Ville Distribution with hanning window")]
        [DataRow("gauss", 51, "gen_kernel_swvd_gauss51_time", DisplayName = "Smoothed Wigner-Ville Distribution with gauss window")]
        public void TestSWVDKernelInTimeDomain(string windowType, int windowLength, string valueList)
        {
            var window = Tools.CreateWindow(windowType);
            Assert.IsNotNull(window);

            var resources = Resources.ResourceManager;
            var mat = new MatReader((byte[])resources.GetObject(valueList));
            var expectedValues = CreateMatrix.DenseOfArray(mat.Read<double[,]>("g"));

            var kernel = new SmoothedWignerVilleDistribution
            {
                Window = window(windowLength),
                UseDopplerDomain = false
            };

            var actualValues = CreateMatrix.DenseOfArray(kernel.Evaluate(expectedValues.RowCount));

            Assert.AreEqual(expectedValues.RowCount, actualValues.RowCount);
            Assert.AreEqual(expectedValues.ColumnCount, actualValues.ColumnCount);

            Assert.AreEqual(0, (expectedValues - actualValues).RowSums().Sum(), 1e-5);
        }

        [TestMethod]
        [DataRow("bart", 51, "gen_kernel_swvd_bart51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with bartlett window")]
        [DataRow("hann", 51, "gen_kernel_swvd_hann51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with hanning window")]
        [DataRow("gauss", 51, "gen_kernel_swvd_gauss51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with gauss window")]
        public void TestSWVDKernelInDopplerDomain(string windowType, int windowLength, string valueList)
        {
            var window = Tools.CreateWindow(windowType);
            Assert.IsNotNull(window);

            var resources = Resources.ResourceManager;
            var mat = new MatReader((byte[])resources.GetObject(valueList));
            var expectedValues = CreateMatrix.DenseOfArray(mat.Read<double[,]>("g"));

            var kernel = new SmoothedWignerVilleDistribution
            {
                Window = window(windowLength),
                UseDopplerDomain = true
            };

            var actualValues = CreateMatrix.DenseOfArray(kernel.Evaluate(expectedValues.RowCount));

            Assert.AreEqual(expectedValues.RowCount, actualValues.RowCount);
            Assert.AreEqual(expectedValues.ColumnCount, actualValues.ColumnCount);

            Assert.AreEqual(0, (expectedValues - actualValues).RowSums().Sum(), 1e-5);
        }
    }
}
