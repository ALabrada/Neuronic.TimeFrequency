using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Resources;
using Accord.Math.Wavelets;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.TimeFrequency.Testing.Properties;
using Neuronic.TimeFrequency.Wavelets;

namespace Neuronic.TimeFrequency.Testing
{
    [TestClass]
    public class WaveletTests
    {
        [TestMethod]
        [DataRow("morl", "wavefun_morl", DisplayName = "Morlet")]
        [DataRow("mexh", "wavefun_mexh", DisplayName = "Mexican Hat")]
        [DataRow("gaus1", "wavefun_gaus1", DisplayName = "Gaussian order 1")]
        [DataRow("gaus2", "wavefun_gaus2", DisplayName = "Gaussian order 2")]
        [DataRow("gaus3", "wavefun_gaus3", DisplayName = "Gaussian order 3")]
        [DataRow("gaus4", "wavefun_gaus4", DisplayName = "Gaussian order 4")]
        [DataRow("gaus5", "wavefun_gaus5", DisplayName = "Gaussian order 5")]

        public void TestEvaluatingContinuousWavelets(string wavName, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            var expectedValues = Tools.ReadNumbersFrom(valueList).Select(c => (Complex) c).ToList();

            var wavelet = (ContinuousWavelet)Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var actualValues = wavelet.Evaluate(wavelet.Minimum, wavelet.Maximum, expectedValues.Count);

            Tools.AssertAreEqual(expectedValues, actualValues, 1e-7);
        }

        [TestMethod]
        [DataRow("haar", "wavefun_haar", DisplayName = "Haar")]
        [DataRow("db5", "wavefun_db5", DisplayName = "Daubechies order 5")]
        [DataRow("db10", "wavefun_db10", DisplayName = "Daubechies order 10")]
        [DataRow("db20", "wavefun_db20", DisplayName = "Daubechies order 20")]
        [DataRow("bior3.3", "wavefun_bior3_3", DisplayName = "Biorthogonal 3.3")]
        [DataRow("bior2.6", "wavefun_bior2_6", DisplayName = "Biorthogonal 3.6")]
        [DataRow("rbio3.3", "wavefun_rbio3_3", DisplayName = "Reverse biorthogonal 3.3")]
        [DataRow("rbio1.5", "wavefun_rbio1_5", DisplayName = "Reverse biorthogonal 1.5")]
        public void TestEvaluatingOrthogonalWavelets(string wavName, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            List<float> x;
            List<Complex> expectedValues;
            using (var reader = new StringReader(valueList))
            {
                x = Tools.ReadNumbersFrom(reader.ReadLine()).ToList();
                expectedValues = Tools.ReadNumbersFrom(reader.ReadLine()).Select(c => (Complex)c).ToList();
            }

            var wavelet = (OrthogonalWavelet)Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            //var actualValues = new Complex[expectedValues.Count];
            //var signal = wavelet.Evaluate();
            //Array.Copy(signal.Samples, actualValues, signal.Count);
            var actualValues = wavelet.Evaluate();

            Tools.AssertAreEqual(expectedValues, actualValues, 1e-5);
        }

        [TestMethod]
        [DataRow("haar", "wfilters_haar", DisplayName = "Haar")]
        [DataRow("db5", "wfilters_db5", DisplayName = "Daubechies order 5")]
        [DataRow("db10", "wfilters_db10", DisplayName = "Daubechies order 10")]
        [DataRow("db20", "wfilters_db20", DisplayName = "Daubechies order 20")]
        [DataRow("bior3.3", "wfilters_bior3_3", DisplayName = "Biorthogonal 3.3")]
        [DataRow("bior2.6", "wfilters_bior2_6", DisplayName = "Biorthogonal 3.6")]
        [DataRow("rbio3.3", "wfilters_rbio3_3", DisplayName = "Reverse biorthogonal 3.3")]
        [DataRow("rbio1.5", "wfilters_rbio1_5", DisplayName = "Reverse biorthogonal 1.5")]
        public void TestOrthogonalFilters(string wavName, string valueList)
        {
            var wavelet = (OrthogonalWavelet)Wavelets.Wavelets.FromName(wavName);
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            List<double> lowRec, highRec, lowRec2 = null, highRec2 = null;
            using (var reader = new StringReader(valueList))
            {
                lowRec = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double) x).ToList();
                highRec = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double) x).ToList();
                if (wavelet is BiorthogonalWavelet)
                {
                    lowRec2 = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double)x).ToList();
                    highRec2 = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double)x).ToList();
                }
            }
            Tools.AssertAreEqual(lowRec, wavelet.LowReconstructionFilter, 1e-5);
            Tools.AssertAreEqual(highRec, wavelet.HighReconstructionFilter, 1e-5);
            if (wavelet is BiorthogonalWavelet bior)
            {
                Tools.AssertAreEqual(lowRec2, bior.Other.LowReconstructionFilter, 1e-5);
                Tools.AssertAreEqual(highRec2, bior.Other.HighReconstructionFilter, 1e-5);
            }
        }

        [TestMethod]
        [DataRow("morl", 0.8125, DisplayName = "Morlet")]
        [DataRow("mexh", 0.25, DisplayName = "Mexican Hat")]
        [DataRow("gaus1", 0.2, DisplayName = "Gaussian order 1")]
        [DataRow("gaus2", 0.3, DisplayName = "Gaussian order 2")]
        [DataRow("gaus3", 0.4, DisplayName = "Gaussian order 3")]
        [DataRow("gaus4", 0.5, DisplayName = "Gaussian order 4")]
        [DataRow("gaus5", 0.5, DisplayName = "Gaussian order 5")]

        public void TestCentralFrequencyOfContinuousWavelets(string wavName, double expectedValue)
        {
            var wavelet = Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var actualValue = wavelet.CentralFrequency;
            Assert.AreEqual(expectedValue, actualValue, 0.1);
        }

        [TestMethod]
        [DataRow("haar", 0.9990, DisplayName = "Haar")]
        [DataRow("db5", 0.6667, DisplayName = "Daubechies order 5")]
        [DataRow("db10", 0.6842, DisplayName = "Daubechies order 10")]
        [DataRow("db20", 0.6667, DisplayName = "Daubechies order 20")]
        [DataRow("bior3.3", 1d, DisplayName = "Biorthogonal 3.3")]
        [DataRow("rbio3.3", 0.4286, DisplayName = "Reverse biorthogonal 3.3")]
        public void TestCentralFrequencyOfOrthogonalWavelets(string wavName, double expectedValue)
        {
            var wavelet = Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var actualValue = wavelet.CentralFrequency;
            Assert.AreEqual(expectedValue, actualValue, 0.1);
        }
    }
}
