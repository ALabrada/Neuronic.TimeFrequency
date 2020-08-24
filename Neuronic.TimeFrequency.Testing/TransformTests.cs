﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using MathNet.Numerics;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Neuronic.TimeFrequency.Kernels;
using Neuronic.TimeFrequency.Testing.Properties;
using Neuronic.TimeFrequency.Transforms;
using Neuronic.TimeFrequency.Wavelets;

namespace Neuronic.TimeFrequency.Testing
{
    [TestClass]
    public class TransformTests
    {
        [TestMethod]
        [DataRow("gauss", 11, 0, "spectrogram_gauss11_noover", DisplayName = "Gabor transform length=11, no overlap")]
        [DataRow("gauss", 11, 10, "spectrogram_gauss11_over10", DisplayName = "Gabor transform length=11, stride 1")]
        [DataRow("hanning", 255, 0, "spectrogram_hann255_noover", DisplayName = "STFT with hanning window length=255, no overlap")]
        [DataRow("bartlett", 125, 62, "spectrogram_bart125_over62", DisplayName = "STFT with bartlett window length=125, 50% overlap")]
        public void TestShortTimeFourierTransform(string winName, int winLength, int overlap, string valueList)
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
                    expectedValues.Add(Tools.ReadComplexNumbersFrom(line).ToArray());
            }

            var window = Tools.CreateWindow(winName).Invoke(winLength);
            var spectrogram = Spectrogram.Estimate(new Signal<float>(samples), window, overlap);

            Assert.AreEqual(expectedValues.Count, spectrogram.FrequencyCount);
            var actualValues = new List<Complex>(samples.Length);
            for (int i = 0; i < expectedValues.Count; i++)
            {
                actualValues.Clear();
                actualValues.AddRange(spectrogram.EnumerateValuesOfFrequencyAt(i));
                Tools.AssertAreEqual(expectedValues[i], actualValues, 1e-3);
            }
        }

        [TestMethod]
        [DataRow("morl", "mycwt_morl", DisplayName = "Morlet")]
        [DataRow("mexh", "mycwt_mexh", DisplayName = "Mexican Hat")]
        [DataRow("gaus1", "mycwt_gaus1", DisplayName = "Gaussian order 1")]
        [DataRow("gaus5", "mycwt_gaus5", DisplayName = "Gaussian order 5")]

        public void TestContinuousWaveletTransformUsingFFT(string wavName, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] t;
            float[] samples;
            var expectedValues = new List<Complex[]>();
            using (var reader = new StringReader(valueList))
            {
                t = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadComplexNumbersFrom(line).ToArray());
            }
            var scales = Enumerable.Range(1, expectedValues.Count).Select(x => (double)x).ToArray();

            var wavelet = Wavelets.Wavelets.FromName(wavName);
            Assert.IsNotNull(wavelet);
            var cwt = ContinuousWaveletTransform.EstimateUsingFFT(new Signal<float>(samples, fs: 1d/(t[1] - t[0])), wavelet, scales);

            Assert.AreEqual(scales.Length, cwt.Scales.Count);
            var actualValues = new List<Complex>(samples.Length);
            for (int i = 0; i < scales.Length; i++)
            {
                actualValues.Clear();
                actualValues.AddRange(cwt.EnumerateValuesOfScaleAt(i));
                Tools.AssertAreEqual(expectedValues[i].Take(expectedValues[i].Length - 1).ToList(), actualValues, 1e-3);
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
                actualValues.AddRange(cwt.EnumerateValuesOfScaleAt(i));
                Tools.AssertAreEqual(expectedValues[i], actualValues, 1e-3);
            }
        }

        [TestMethod]
        [DataRow("db5", "dwt_db5", DisplayName = "Daubechies order 5")]
        [DataRow("db20", "dwt_db20", DisplayName = "Daubechies order 20")]
        [DataRow("bior3.3", "dwt_bior3_3", DisplayName = "Biorthogonal 3.3")]
        [DataRow("bior2.6", "dwt_bior2_6", DisplayName = "Biorthogonal 2.6")]
        [DataRow("rbio3.3", "dwt_rbio3_3", DisplayName = "Reverse biorthogonal 3.3")]
        [DataRow("rbio1.5", "dwt_rbio1_5", DisplayName = "Reverse biorthogonal 1.5")]

        public void TestDiscreteWaveletTransformWithSymetricPadding(string wavName, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            double[] samples;
            double[] a, d;
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double)x).ToArray();
                a = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double) x).ToArray();
                d = Tools.ReadNumbersFrom(reader.ReadLine()).Select(x => (double)x).ToArray();
            }

            var wavelet = Wavelets.Wavelets.FromName(wavName) as OrthogonalWavelet;
            Assert.IsNotNull(wavelet);
            var dwt = DiscreteWaveletTransform.Estimate(new Signal<double>(samples), wavelet, null);

            Tools.AssertAreEqual(a, dwt.Approximation, 1e-3);
            Tools.AssertAreEqual(d, dwt.Detail, 1e-3);
        }

        [TestMethod]
        [DataRow(30d, "dtfd_cw_30", DisplayName = "Choi-Williams sigma=30")]
        [DataRow(100d, "dtfd_cw_100", DisplayName = "Choi-Williams sigma=100")]
        public void TestChoiWilliamsTimeFrequencyDistribution(double sigma, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] samples;
            var expectedValues = new List<double[]>();
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(x => (double)x).ToArray());
            }

            var kernel = new ChoiWilliamsDistribution {Sigma = sigma};
            var tfd = TimeFrequencyDistribution.Estimate(new Signal<float>(samples, fs: 10d), kernel);

            Assert.AreEqual(expectedValues.Count, tfd.SampleCount);
            Assert.AreEqual(tfd.FrequencyCount, samples.Length);
            var actualValues = new List<double>(samples.Length);
            for (int i = 0; i < expectedValues.Count; i++)
            {
                actualValues.Clear();
                actualValues.AddRange(tfd.EnumerateValuesOfTimeAt(i));
                Tools.AssertAreEqual(expectedValues[i], actualValues, 1e-5);
            }
        }

        [TestMethod]
        [DataRow("ParabEmd_128", DisplayName = "EMD 128 samples with default parameters.")]
        public void TestEmpiricalModeDecomposition(string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] samples;
            var expectedValues = new List<double[]>();
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(x => (double)x).ToArray());
            }

            var emd = EmpiricalModeDecomposition.Estimate(new Signal<float>(samples));

            Assert.AreEqual(expectedValues.Count, emd.Count);
            for (int i = 0; i < expectedValues.Count; i++)
            {
                Tools.AssertAreEqual(expectedValues[i], emd[i], 1e-4*(1<<i));
            }
        }

        [TestMethod]
        [DataRow(10d, "HHT_128", DisplayName = "HHT 128 samples with default parameters.")]
        public void TestHilbertHuangTransform(double fs, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] samples;
            var expectedValues = new List<double[]>();
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
                string line;
                while ((line = reader.ReadLine()) != null)
                    expectedValues.Add(Tools.ReadNumbersFrom(line).Select(x => (double)x).Take(samples.Length - 1).ToArray());
            }

            var emd = EmpiricalModeDecomposition.Estimate(new Signal<float>(samples, fs: fs));
            var hht = emd.HilbertSpectralAnalysis();

            Assert.AreEqual(expectedValues.Count, hht.Count);
            for (int i = 0; i < expectedValues.Count; i++)
            {
                Tools.AssertAreEqual(expectedValues[i], hht[i].Frequency, 1e-4 * (2 << i));
            }
        }

        [TestMethod]
        [DataRow(10d, "HHT_128", DisplayName = "HHT 128 samples with default parameters.")]
        public void TestHilbertHuangTransformConversionToSpectrogram(double fs, string valueList)
        {
            var resources = Resources.ResourceManager;
            valueList = resources.GetString(valueList) ?? valueList;
            float[] samples;
            using (var reader = new StringReader(valueList))
            {
                samples = Tools.ReadNumbersFrom(reader.ReadLine()).ToArray();
            }

            var emd = EmpiricalModeDecomposition.Estimate(new Signal<float>(samples, fs: fs));
            var hht = emd.HilbertSpectralAnalysis();
            var spectrogram = hht.GetSpectrogram();

            for (int i = 0; i < hht.SampleCount; i++)
            {
                var randomValues = new List<int>(10);
                randomValues.AddRange(Generate.UniformSequence().Take(10).Select(j => (int)(j * (spectrogram.FrequencyCount - 1))).Distinct());

                foreach (var component in hht)
                {
                    if (i - component.FrequencyOffset >= 0 && i - component.FrequencyOffset < component.Frequency.Count)
                    {
                        var frequency = component.Frequency[i - component.FrequencyOffset];
                        var freqIndex = spectrogram.FindClosestIndexOfFrequency(frequency);
                        randomValues.Remove(freqIndex);
                        randomValues.Remove(freqIndex - 1);
                        randomValues.Remove(freqIndex + 1);

                        var expectedValue = component.Amplitude[i];
                        var actualValue = spectrogram[i, freqIndex];
                        Assert.IsTrue(actualValue >= expectedValue);
                    }
                }

                foreach (var freqIndex in randomValues)
                    Assert.AreEqual(0, spectrogram[i, freqIndex]);
            }
        }
    }
}
