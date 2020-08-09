using System;
using Accord.Math.Wavelets;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace Neuronic.TimeFrequency.Testing
{
    [TestClass]
    public class WaveletTests
    {
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
