using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.TimeFrequency.Testing.Properties;

namespace Neuronic.TimeFrequency.Testing
{
    [TestClass]
    public class TransformTests
    {
        [TestMethod]
        [DataRow("morl", "cwt_morl", DisplayName = "Morlet")]
        [DataRow("mexh", "cwt_mexh", DisplayName = "Mexican Hat")]
        [DataRow("gaus1", "cwt_gaus1", DisplayName = "Gaussian order 1")]
        [DataRow("gaus5", "cwt_gaus5", DisplayName = "Gaussian order 5")]

        public void TestContinuousWaveletTransformUsingFFT(string wavName, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] samples;
            var expectedValues = new List<Complex[]>();
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(x => (Complex)x).ToArray());
            }
            var scales = Enumerable.Range(1, expectedValues.Count).Select(x => (double)x).ToArray();

            var wavelet = Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var cwt = ContinuousWaveletTransform.EstimateUsingFFT(new Signal<float>(samples), wavelet, scales);

            Assert.AreEqual(scales.Length, cwt.Scales.Count);
            var actualValues = new List<Complex>(samples.Length);
            for (int i = 0; i < scales.Length; i++)
            {
                actualValues.Clear();
                actualValues.AddRange(cwt.EnumerateScale(i));
                Tools.AssertAreEqual(expectedValues[i], actualValues, 1e-3);
            }            
        }

        [TestMethod]
        [DataRow("morl", "cwt_morl", DisplayName = "Morlet")]
        [DataRow("mexh", "cwt_mexh", DisplayName = "Mexican Hat")]
        [DataRow("gaus1", "cwt_gaus1", DisplayName = "Gaussian order 1")]
        [DataRow("gaus5", "cwt_gaus5", DisplayName = "Gaussian order 5")]
        [DataRow("db5", "cwt_db5", DisplayName = "Daubechies order 5")]
        [DataRow("db20", "cwt_db20", DisplayName = "Daubechies order 20")]
        [DataRow("bior3.3", "cwt_bior3_3", DisplayName = "Biorthogonal 3.3")]
        [DataRow("rbio3.3", "cwt_rbio3_3", DisplayName = "Reverse biorthogonal 3.3")]

        public void TestContinuousWaveletTransformUsingConvolution(string wavName, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] samples;
            var expectedValues = new List<Complex[]>();
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(x => (Complex)x).ToArray());
            }
            var scales = Enumerable.Range(1, expectedValues.Count).Select(x => (double)x).ToArray();

            var wavelet = Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var cwt = ContinuousWaveletTransform.EstimateUsingConvolutions(new Signal<float>(samples), wavelet, scales);

            Assert.AreEqual(scales.Length, cwt.Scales.Count);
            var actualValues = new List<Complex>(samples.Length);
            for (int i = 0; i < scales.Length; i++)
            {
                actualValues.Clear();
                actualValues.AddRange(cwt.EnumerateScale(i));
                Tools.AssertAreEqual(expectedValues[i], actualValues, 1e-3);
            }
        }
    }
}
