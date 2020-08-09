using System;
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
            var resources = new ResourceManager("Neuronic.TimeFrequency.Testing.Properties.Resources",
                typeof(Resources).Assembly);
            valueList = resources.GetString(valueList) ?? valueList;
            var expectedValues = Tools.ReadNumbersFrom(valueList).Select(c => (Complex) c).ToList();

            var wavelet = (ContinuousWavelet)Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var actualValues = wavelet.Evaluate(wavelet.Minimum, wavelet.Maximum, expectedValues.Count);

            Tools.AssertAreEqual(expectedValues, actualValues, 1e-7);
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
