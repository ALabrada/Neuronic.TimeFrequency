﻿using System;
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

        [TestMethod]
        [DataRow("delta", 5, "0	0	1	0	0", DisplayName = "Delta window")]
        [DataRow("rect", 5, "1	1	1	1	1", DisplayName = "Rectangular window")]
        [DataRow("bart", 5, "0    0.5000    1.0000    0.5000         0", DisplayName = "Bartlett window")]
        [DataRow("hamm", 5, "0.0800    0.5400    1.0000    0.5400    0.0800", DisplayName = "Hamming window")]
        [DataRow("hann", 5, "0.2500    0.7500    1.0000    0.7500    0.2500", DisplayName = "Hanning window")]
        [DataRow("gauss", 5, "0.0439    0.4578    1.0000    0.4578    0.0439", DisplayName = "Gauss window")]
        public void TestWindowFunctions(string windowType, int windowLength, string valueList)
        {
            var window = Tools.CreateWindow(windowType);
            Assert.IsNotNull(window);

            var expectedValues = Tools.ReadNumbersFrom(valueList).Select(c => (double) c).ToArray();

            var actualValues = window.Evaluate(windowLength);

            Tools.AssertAreEqual(expectedValues, actualValues, 1e-4);
        }

        [TestMethod]
        [DataRow("delta", 51, "gen_kernel_swvd_delta51_time", DisplayName = "Smoothed Wigner-Ville Distribution with delta window")]
        [DataRow("rect", 51, "gen_kernel_swvd_rect51_time", DisplayName = "Smoothed Wigner-Ville Distribution with rectangular window")]
        [DataRow("bart", 51, "gen_kernel_swvd_bart51_time", DisplayName = "Smoothed Wigner-Ville Distribution with bartlett window")]
        [DataRow("hamm", 51, "gen_kernel_swvd_hamm51_time", DisplayName = "Smoothed Wigner-Ville Distribution with hamming window")]
        [DataRow("hann", 51, "gen_kernel_swvd_hann51_time", DisplayName = "Smoothed Wigner-Ville Distribution with hanning window")]
        [DataRow("gauss", 51, "gen_kernel_swvd_gauss51_time", DisplayName = "Smoothed Wigner-Ville Distribution with gauss window")]
        public void TestSWVDKernelInTimeDomain(string windowType, int windowLength, string valueList)
        {
            var window = Tools.CreateWindow(windowType);
            Assert.IsNotNull(window);

            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            var expectedValues = new List<double[]>();
            using (var reader = new StringReader(valueList))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(c => (double)c).ToArray());
            }

            var kernel = new SmoothedWignerVilleDistribution
            {
                Window = window.Evaluate(windowLength),
                UseDopplerDomain = false
            };

            var actualValues = kernel.Evaluate(expectedValues.Count);

            Assert.AreEqual(expectedValues.Count, actualValues.GetLength(0));
            for (int i = 0; i < expectedValues.Count; i++)
            {
                var row = new double[actualValues.GetLength(1)];
                Buffer.BlockCopy(actualValues, i * row.Length * sizeof(double), row, 0, row.Length * sizeof(double));
                Tools.AssertAreEqual(expectedValues[i], row, 1e-3);
            }
        }

        [TestMethod]
        [DataRow("delta", 51, "gen_kernel_swvd_delta51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with delta window")]
        [DataRow("rect", 51, "gen_kernel_swvd_rect51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with rectangular window")]
        [DataRow("bart", 51, "gen_kernel_swvd_bart51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with bartlett window")]
        [DataRow("hamm", 51, "gen_kernel_swvd_hamm51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with hamming window")]
        [DataRow("hann", 51, "gen_kernel_swvd_hann51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with hanning window")]
        [DataRow("gauss", 51, "gen_kernel_swvd_gauss51_doppler", DisplayName = "Smoothed Wigner-Ville Distribution with gauss window")]
        public void TestSWVDKernelInDopplerDomain(string windowType, int windowLength, string valueList)
        {
            var window = Tools.CreateWindow(windowType);
            Assert.IsNotNull(window);

            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            var expectedValues = new List<double[]>();
            using (var reader = new StringReader(valueList))
            {
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(c => (double)c).ToArray());
            }

            var kernel = new SmoothedWignerVilleDistribution
            {
                Window = window.Evaluate(windowLength),
                UseDopplerDomain = true
            };

            var actualValues = kernel.Evaluate(expectedValues.Count);

            Assert.AreEqual(expectedValues.Count, actualValues.GetLength(0));
            for (int i = 0; i < expectedValues.Count; i++)
            {
                var row = new double[actualValues.GetLength(1)];
                Buffer.BlockCopy(actualValues, i * row.Length * sizeof(double), row, 0, row.Length * sizeof(double));
                Tools.AssertAreEqual(expectedValues[i], row, 1e-3);
            }
        }
    }
}