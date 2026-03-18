const palette = {
  sc: "#0f766e",
  mrc: "#b45309",
  sskmSc: "#7c3aed",
  sskmMrc: "#b42318",
  siso: "#0f766e",
};

const dom = {
  scenario: document.getElementById("scenario"),
  fftSize: document.getElementById("fftSize"),
  giRatio: document.getElementById("giRatio"),
  modulationOrder: document.getElementById("modulationOrder"),
  multiPath: document.getElementById("multiPath"),
  iterations: document.getElementById("iterations"),
  snrStart: document.getElementById("snrStart"),
  snrEnd: document.getElementById("snrEnd"),
  snrStep: document.getElementById("snrStep"),
  seed: document.getElementById("seed"),
  runButton: document.getElementById("runButton"),
  fillPresetButton: document.getElementById("fillPresetButton"),
  chartCanvas: document.getElementById("chartCanvas"),
  legend: document.getElementById("legend"),
  summaryCards: document.getElementById("summaryCards"),
  resultsTable: document.getElementById("resultsTable"),
  resultsHeader: document.getElementById("resultsHeader"),
  statusText: document.getElementById("statusText"),
};

class RNG {
  constructor(seed) {
    this.state = (seed >>> 0) || 1;
  }

  next() {
    this.state = (1664525 * this.state + 1013904223) >>> 0;
    return this.state / 4294967296;
  }

  normalPair() {
    const u1 = Math.max(this.next(), 1e-12);
    const u2 = this.next();
    const mag = Math.sqrt(-2 * Math.log(u1));
    const angle = 2 * Math.PI * u2;
    return [mag * Math.cos(angle), mag * Math.sin(angle)];
  }
}

function complex(re = 0, im = 0) {
  return { re, im };
}

function add(a, b) {
  return { re: a.re + b.re, im: a.im + b.im };
}

function sub(a, b) {
  return { re: a.re - b.re, im: a.im - b.im };
}

function mul(a, b) {
  return {
    re: a.re * b.re - a.im * b.im,
    im: a.re * b.im + a.im * b.re,
  };
}

function scale(a, value) {
  return { re: a.re * value, im: a.im * value };
}

function conj(a) {
  return { re: a.re, im: -a.im };
}

function abs2(a) {
  return a.re * a.re + a.im * a.im;
}

function div(a, b) {
  const denom = abs2(b) + 1e-12;
  return {
    re: (a.re * b.re + a.im * b.im) / denom,
    im: (a.im * b.re - a.re * b.im) / denom,
  };
}

function fft(input) {
  const n = input.length;
  if (n <= 1) {
    return input.slice();
  }

  const even = [];
  const odd = [];
  for (let i = 0; i < n; i += 1) {
    (i % 2 === 0 ? even : odd).push(input[i]);
  }

  const fftEven = fft(even);
  const fftOdd = fft(odd);
  const result = new Array(n);

  for (let k = 0; k < n / 2; k += 1) {
    const angle = (-2 * Math.PI * k) / n;
    const twiddle = complex(Math.cos(angle), Math.sin(angle));
    const t = mul(twiddle, fftOdd[k]);
    result[k] = add(fftEven[k], t);
    result[k + n / 2] = sub(fftEven[k], t);
  }

  return result;
}

function ifft(input) {
  const conjugated = input.map(conj);
  const transformed = fft(conjugated).map(conj);
  return transformed.map((value) => scale(value, 1 / input.length));
}

function convolve(signal, channel) {
  const out = Array.from(
    { length: signal.length + channel.length - 1 },
    () => complex(0, 0),
  );

  for (let i = 0; i < signal.length; i += 1) {
    for (let j = 0; j < channel.length; j += 1) {
      out[i + j] = add(out[i + j], mul(signal[i], channel[j]));
    }
  }

  return out;
}

function createSnrRange(start, end, step) {
  const values = [];
  for (let snr = start; snr <= end + 1e-9; snr += step) {
    values.push(Number(snr.toFixed(6)));
  }
  return values;
}

function randomBits(length, rng) {
  const bits = new Array(length);
  for (let i = 0; i < length; i += 1) {
    bits[i] = rng.next() > 0.5 ? 1 : 0;
  }
  return bits;
}

function rayleighChannel(multiPath, rng) {
  const profile = [];
  let total = 0;
  for (let i = 0; i < multiPath; i += 1) {
    const power = Math.exp(-(i + 1) / 5);
    profile.push(power);
    total += power;
  }

  return profile.map((power) => {
    const [n1, n2] = rng.normalPair();
    const sigma = Math.sqrt(0.5 * (power / total));
    return complex(n1 * sigma, n2 * sigma);
  });
}

function addAwgn(signal, snrDb, rng) {
  const noisePower = 10 ** (-snrDb / 10);
  const sigma = Math.sqrt(noisePower / 2);
  return signal.map((sample) => {
    const [n1, n2] = rng.normalPair();
    return add(sample, complex(n1 * sigma, n2 * sigma));
  });
}

function qpskMod(bits) {
  const symbols = [];
  for (let i = 0; i < bits.length; i += 2) {
    const imagBit = bits[i];
    const realBit = bits[i + 1];
    const imag = imagBit ? 1 : -1;
    const real = realBit ? 1 : -1;
    symbols.push(complex(real * 0.7071, imag * 0.7071));
  }
  return symbols;
}

function qpskDemod(symbols) {
  const bits = [];
  for (const symbol of symbols) {
    bits.push(symbol.im > 0 ? 1 : 0);
    bits.push(symbol.re > 0 ? 1 : 0);
  }
  return bits;
}

function qam16Mod(bits) {
  const symbols = [];
  for (let i = 0; i < bits.length; i += 4) {
    const b1 = bits[i];
    const b2 = bits[i + 1];
    const b3 = bits[i + 2];
    const b4 = bits[i + 3];
    const first = b1 * 4 - 2;
    const second = (b2 !== b1 ? 1 : 0) * 2 - 1;
    const third = b3 * 4 - 2;
    const fourth = (b4 !== b3 ? 1 : 0) * 2 - 1;
    symbols.push(complex((first + second) * 0.3162, (third + fourth) * 0.3162));
  }
  return symbols;
}

function qam16Demod(symbols) {
  const bits = [];
  for (const symbol of symbols) {
    bits.push(symbol.re > 0 ? 1 : 0);
    bits.push(Math.abs(symbol.re) < 0.6325 ? 1 : 0);
    bits.push(symbol.im > 0 ? 1 : 0);
    bits.push(Math.abs(symbol.im) < 0.6325 ? 1 : 0);
  }
  return bits;
}

function baseMod(bits, modulationOrder) {
  return modulationOrder === 2 ? qpskMod(bits) : qam16Mod(bits);
}

function baseDemod(symbols, modulationOrder) {
  return modulationOrder === 2 ? qpskDemod(symbols) : qam16Demod(symbols);
}

function addCyclicPrefix(ofdm, giSize) {
  return ofdm.slice(ofdm.length - giSize).concat(ofdm);
}

function removeCyclicPrefix(rx, giSize, fftSize) {
  return rx.slice(giSize, giSize + fftSize);
}

function normalizeIfft(symbols) {
  const factor = Math.sqrt(symbols.length);
  return ifft(symbols).map((value) => scale(value, factor));
}

function normalizeFft(samples) {
  const factor = Math.sqrt(samples.length);
  return fft(samples).map((value) => scale(value, 1 / factor));
}

function channelResponse(channel, fftSize) {
  const padded = Array.from({ length: fftSize }, (_, index) => channel[index] || complex(0, 0));
  return fft(padded);
}

function countBitErrors(a, b) {
  let count = 0;
  for (let i = 0; i < a.length; i += 1) {
    if (a[i] !== b[i]) {
      count += 1;
    }
  }
  return count;
}

function equalize(received, channel) {
  return received.map((value, index) => div(value, channel[index]));
}

function selectCombine(receivedA, receivedB, channelA, channelB) {
  const y = [];
  const h = [];
  for (let i = 0; i < receivedA.length; i += 1) {
    if (abs2(channelA[i]) >= abs2(channelB[i])) {
      y.push(receivedA[i]);
      h.push(channelA[i]);
    } else {
      y.push(receivedB[i]);
      h.push(channelB[i]);
    }
  }
  return { y, h };
}

function mrcCombine(receivedA, receivedB, channelA, channelB) {
  const combined = [];
  const denom = [];
  for (let i = 0; i < receivedA.length; i += 1) {
    const weighted = add(
      mul(receivedA[i], conj(channelA[i])),
      mul(receivedB[i], conj(channelB[i])),
    );
    combined.push(weighted);
    denom.push(complex(abs2(channelA[i]) + abs2(channelB[i]), 0));
  }
  return { y: combined, h: denom };
}

function encodeSskm(bits) {
  const encoded = new Array(bits.length * 2);
  for (let i = 0; i < bits.length; i += 1) {
    if (bits[i] === 1) {
      encoded[2 * i] = 1;
      encoded[2 * i + 1] = 0;
    } else {
      encoded[2 * i] = 0;
      encoded[2 * i + 1] = 1;
    }
  }
  return encoded;
}

function decodeSskm(bits) {
  const decoded = new Array(bits.length / 2);
  for (let i = 0; i < decoded.length; i += 1) {
    decoded[i] = bits[2 * i] > bits[2 * i + 1] ? 1 : 0;
  }
  return decoded;
}

function simulateAssignment2(config, rng) {
  const fftSize = config.fftSize;
  const giSize = Math.floor(fftSize * config.giRatio);
  const dataSize = fftSize * config.modulationOrder;
  const snrValues = createSnrRange(config.snrStart, config.snrEnd, config.snrStep);
  const ber = [];

  for (const snr of snrValues) {
    let errors = 0;
    for (let iter = 0; iter < config.iterations; iter += 1) {
      const data = randomBits(dataSize, rng);
      const modulated = baseMod(data, config.modulationOrder);
      const ofdm = normalizeIfft(modulated);
      const withCp = addCyclicPrefix(ofdm, giSize);
      const channel = rayleighChannel(config.multiPath, rng);
      const rx = addAwgn(convolve(withCp, channel), snr, rng);
      const trimmed = removeCyclicPrefix(rx, giSize, fftSize);
      const freq = normalizeFft(trimmed);
      const h = channelResponse(channel, fftSize);
      const equalized = equalize(freq, h);
      const demod = baseDemod(equalized, config.modulationOrder).slice(0, data.length);
      errors += countBitErrors(data, demod);
    }
    ber.push(errors / (dataSize * config.iterations));
  }

  return {
    snrValues,
    series: [
      {
        key: "siso",
        label: "SISO-OFDM",
        color: palette.siso,
        values: ber,
      },
    ],
  };
}

function simulateAssignment3(config, rng) {
  const fftSize = config.fftSize;
  const fftSizeSskm = fftSize * 2;
  const giSize = Math.floor(fftSize * config.giRatio);
  const giSizeSskm = Math.floor(fftSizeSskm * config.giRatio);
  const dataSize = fftSize * config.modulationOrder;
  const dataSizeSskm = dataSize * 2;
  const snrValues = createSnrRange(config.snrStart, config.snrEnd, config.snrStep);
  const berSc = [];
  const berMrc = [];
  const berSskmSc = [];
  const berSskmMrc = [];

  for (const snr of snrValues) {
    let errorsSc = 0;
    let errorsMrc = 0;
    let errorsSskmSc = 0;
    let errorsSskmMrc = 0;

    for (let iter = 0; iter < config.iterations; iter += 1) {
      const data = randomBits(dataSize, rng);
      const modulated = baseMod(data, config.modulationOrder);
      const ofdm = normalizeIfft(modulated);
      const tx = addCyclicPrefix(ofdm, giSize);

      const h1 = rayleighChannel(config.multiPath, rng);
      const h2 = rayleighChannel(config.multiPath, rng);
      const y1 = addAwgn(convolve(tx, h1), snr, rng);
      const y2 = addAwgn(convolve(tx, h2), snr, rng);

      const rx1 = normalizeFft(removeCyclicPrefix(y1, giSize, fftSize));
      const rx2 = normalizeFft(removeCyclicPrefix(y2, giSize, fftSize));
      const H1 = channelResponse(h1, fftSize);
      const H2 = channelResponse(h2, fftSize);

      const sc = selectCombine(rx1, rx2, H1, H2);
      const mrc = mrcCombine(rx1, rx2, H1, H2);
      const demodSc = baseDemod(equalize(sc.y, sc.h), config.modulationOrder).slice(0, data.length);
      const demodMrc = baseDemod(equalize(mrc.y, mrc.h), config.modulationOrder).slice(0, data.length);

      errorsSc += countBitErrors(data, demodSc);
      errorsMrc += countBitErrors(data, demodMrc);

      const dataSskm = encodeSskm(data);
      const modulatedSskm = baseMod(dataSskm, config.modulationOrder);
      const ofdmSskm = normalizeIfft(modulatedSskm);
      const txSskm = addCyclicPrefix(ofdmSskm, giSizeSskm);
      const y1Sskm = addAwgn(convolve(txSskm, h1), snr, rng);
      const y2Sskm = addAwgn(convolve(txSskm, h2), snr, rng);
      const rx1Sskm = normalizeFft(removeCyclicPrefix(y1Sskm, giSizeSskm, fftSizeSskm));
      const rx2Sskm = normalizeFft(removeCyclicPrefix(y2Sskm, giSizeSskm, fftSizeSskm));
      const H1Sskm = channelResponse(h1, fftSizeSskm);
      const H2Sskm = channelResponse(h2, fftSizeSskm);

      const scSskm = selectCombine(rx1Sskm, rx2Sskm, H1Sskm, H2Sskm);
      const mrcSskm = mrcCombine(rx1Sskm, rx2Sskm, H1Sskm, H2Sskm);
      const decodedSc = decodeSskm(
        baseDemod(equalize(scSskm.y, scSskm.h), config.modulationOrder).slice(0, dataSizeSskm),
      );
      const decodedMrc = decodeSskm(
        baseDemod(equalize(mrcSskm.y, mrcSskm.h), config.modulationOrder).slice(0, dataSizeSskm),
      );

      errorsSskmSc += countBitErrors(data, decodedSc);
      errorsSskmMrc += countBitErrors(data, decodedMrc);
    }

    berSc.push(errorsSc / (dataSize * config.iterations));
    berMrc.push(errorsMrc / (dataSize * config.iterations));
    berSskmSc.push(errorsSskmSc / (dataSize * config.iterations));
    berSskmMrc.push(errorsSskmMrc / (dataSize * config.iterations));
  }

  return {
    snrValues,
    series: [
      { key: "sc", label: "SC", color: palette.sc, values: berSc },
      { key: "mrc", label: "MRC", color: palette.mrc, values: berMrc },
      { key: "sskmSc", label: "SSKM-SC", color: palette.sskmSc, values: berSskmSc },
      { key: "sskmMrc", label: "SSKM-MRC", color: palette.sskmMrc, values: berSskmMrc },
    ],
  };
}

function getConfig() {
  const config = {
    scenario: dom.scenario.value,
    fftSize: Number(dom.fftSize.value),
    giRatio: Number(dom.giRatio.value),
    modulationOrder: Number(dom.modulationOrder.value),
    multiPath: Number(dom.multiPath.value),
    iterations: Number(dom.iterations.value),
    snrStart: Number(dom.snrStart.value),
    snrEnd: Number(dom.snrEnd.value),
    snrStep: Number(dom.snrStep.value),
    seed: Number(dom.seed.value),
  };

  if (config.fftSize <= 0 || (config.fftSize & (config.fftSize - 1)) !== 0) {
    throw new Error("FFT Size must be a power of two.");
  }
  if (config.snrEnd < config.snrStart) {
    throw new Error("SNR End must be greater than or equal to SNR Start.");
  }
  if (config.snrStep <= 0) {
    throw new Error("SNR Step must be greater than zero.");
  }
  if (config.modulationOrder !== 2 && config.modulationOrder !== 4) {
    throw new Error("This web simulator currently supports QPSK and 16QAM only.");
  }

  return config;
}

function formatBer(value) {
  return value < 0.001 ? value.toExponential(2) : value.toFixed(4);
}

function renderSummary(result) {
  dom.summaryCards.innerHTML = "";

  const cards = result.series.map((serie) => {
    const best = Math.min(...serie.values);
    const worst = Math.max(...serie.values);
    return `
      <div class="summary-card">
        <div>${serie.label}</div>
        <strong>${formatBer(best)}</strong>
        <div>Best BER</div>
        <div>Worst ${formatBer(worst)}</div>
      </div>
    `;
  }).join("");

  dom.summaryCards.innerHTML = cards;
}

function renderTable(result) {
  dom.resultsHeader.textContent = result.series.map((serie) => serie.label).join(" / ");
  const rows = result.snrValues.map((snr, index) => {
    const values = result.series.map((serie) => formatBer(serie.values[index])).join(" / ");
    return `<tr><td>${snr}</td><td>${values}</td></tr>`;
  }).join("");

  dom.resultsTable.innerHTML = rows;
}

function renderLegend(series) {
  dom.legend.innerHTML = series.map((serie) => `
    <div class="legend-item">
      <span class="legend-swatch" style="background:${serie.color}"></span>
      <span>${serie.label}</span>
    </div>
  `).join("");
}

function drawChart(result) {
  const canvas = dom.chartCanvas;
  const ctx = canvas.getContext("2d");
  const width = canvas.width;
  const height = canvas.height;
  const pad = { left: 72, right: 24, top: 22, bottom: 52 };
  const minBer = 1e-4;
  const maxBer = 1;
  const plotWidth = width - pad.left - pad.right;
  const plotHeight = height - pad.top - pad.bottom;

  ctx.clearRect(0, 0, width, height);
  ctx.fillStyle = "#fffefb";
  ctx.fillRect(0, 0, width, height);

  ctx.strokeStyle = "rgba(31, 26, 22, 0.12)";
  ctx.lineWidth = 1;
  for (let i = 0; i <= 4; i += 1) {
    const y = pad.top + (plotHeight / 4) * i;
    ctx.beginPath();
    ctx.moveTo(pad.left, y);
    ctx.lineTo(width - pad.right, y);
    ctx.stroke();
  }

  const xMin = result.snrValues[0];
  const xMax = result.snrValues[result.snrValues.length - 1];

  const xFor = (value) => {
    if (xMax === xMin) {
      return pad.left + plotWidth / 2;
    }
    return pad.left + ((value - xMin) / (xMax - xMin)) * plotWidth;
  };

  const yFor = (value) => {
    const clamped = Math.max(minBer, Math.min(maxBer, value));
    const logValue = Math.log10(clamped);
    const ratio = (Math.log10(maxBer) - logValue) / (Math.log10(maxBer) - Math.log10(minBer));
    return pad.top + ratio * plotHeight;
  };

  ctx.strokeStyle = "#34291d";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(pad.left, pad.top);
  ctx.lineTo(pad.left, height - pad.bottom);
  ctx.lineTo(width - pad.right, height - pad.bottom);
  ctx.stroke();

  ctx.fillStyle = "#66584b";
  ctx.font = "12px Segoe UI";
  const yTicks = [1, 1e-1, 1e-2, 1e-3, 1e-4];
  for (const tick of yTicks) {
    const y = yFor(tick);
    ctx.fillText(tick.toExponential(0), 16, y + 4);
  }

  for (const snr of result.snrValues) {
    const x = xFor(snr);
    ctx.fillText(String(snr), x - 8, height - 20);
  }

  ctx.save();
  ctx.translate(18, height / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText("BER (log scale)", 0, 0);
  ctx.restore();
  ctx.fillText("SNR (dB)", width / 2 - 20, height - 8);

  for (const serie of result.series) {
    ctx.strokeStyle = serie.color;
    ctx.fillStyle = serie.color;
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    serie.values.forEach((value, index) => {
      const x = xFor(result.snrValues[index]);
      const y = yFor(value);
      if (index === 0) {
        ctx.moveTo(x, y);
      } else {
        ctx.lineTo(x, y);
      }
    });
    ctx.stroke();

    for (let index = 0; index < serie.values.length; index += 1) {
      const x = xFor(result.snrValues[index]);
      const y = yFor(serie.values[index]);
      ctx.beginPath();
      ctx.arc(x, y, 4, 0, Math.PI * 2);
      ctx.fill();
    }
  }
}

async function runSimulation() {
  try {
    dom.runButton.disabled = true;
    dom.statusText.textContent = "Running...";

    await new Promise((resolve) => setTimeout(resolve, 20));

    const config = getConfig();
    const rng = new RNG(config.seed);
    const result = config.scenario === "assignment2"
      ? simulateAssignment2(config, rng)
      : simulateAssignment3(config, rng);

    renderLegend(result.series);
    renderSummary(result);
    renderTable(result);
    drawChart(result);
    dom.statusText.textContent = `Done · ${result.snrValues.length} SNR points`;
  } catch (error) {
    dom.statusText.textContent = error.message;
    dom.legend.innerHTML = "";
    dom.summaryCards.innerHTML = "";
    dom.resultsTable.innerHTML = "";
    const ctx = dom.chartCanvas.getContext("2d");
    ctx.clearRect(0, 0, dom.chartCanvas.width, dom.chartCanvas.height);
  } finally {
    dom.runButton.disabled = false;
  }
}

function applyAssignment3Preset() {
  dom.scenario.value = "assignment3";
  dom.fftSize.value = 128;
  dom.giRatio.value = "0.25";
  dom.modulationOrder.value = 4;
  dom.multiPath.value = 7;
  dom.iterations.value = 300;
  dom.snrStart.value = 0;
  dom.snrEnd.value = 30;
  dom.snrStep.value = 3;
}

dom.runButton.addEventListener("click", runSimulation);
dom.fillPresetButton.addEventListener("click", () => {
  applyAssignment3Preset();
  runSimulation();
});

applyAssignment3Preset();
runSimulation();
