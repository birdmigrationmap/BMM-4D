"use strict";
// Web map viewer of BMM-4D inversions (webmap/, served by serve.py).
// The exported frames (bmm4d.export_webmap) are column densities on a grid
// regular in longitude and Web Mercator y, drawn by a WebGL custom layer.
// Between two 30-min frames the layer morphs as the MATLAB video: frame i
// advected forward by tau * dt and frame i+1 backward by (1 - tau) * dt
// along the birds' mean wind (two-stage, wind at the half-way point),
// blended with weights 1 - tau and tau. The flow overlay (particles or
// arrows) uses the same wind, and particles move at the birds' own speed.

const EARTH = 40075016.686;               // m, circumference (Web Mercator unit)
const params = new URLSearchParams(location.search);
const config = window.WEBMAP || {};       // set by the shared (Artifact) build, empty in dev
const el = (id) => document.getElementById(id);

const state = {
  meta: null, frame: 0, playing: false, speed: 1800, morph: true, night: true, lowfade: false,
  opacity: 0.85, climit: 200, flow: "particles", particleCount: 4000, radars: false,
  cols: [], winds: [], lastNow: null,
};
let map = null;
let lib = null;

// ---------------------------------------------------------------- helpers
const pad = (i) => String(i).padStart(2, "0");
const latOfMercY = (y) => Math.atan(Math.sinh(Math.PI * (1 - 2 * y))) * 180 / Math.PI;
// Mercator units per metre at Mercator y: 1 / (EARTH cos(lat)), with cos(lat) = 1 / cosh(pi (1 - 2 y))
const mercPerMetre = (y) => Math.cosh(Math.PI * (1 - 2 * y)) / EARTH;
const W = new Float64Array(2);            // scratch wind (m/s, east and north)

function loadScript(src) {
  return new Promise((resolve, reject) => {
    const s = document.createElement("script");
    s.src = src; s.onload = resolve; s.onerror = reject;
    document.head.appendChild(s);
  });
}
function loadImage(src) {
  return new Promise((resolve, reject) => {
    const image = new Image();
    image.onload = () => resolve(image);
    image.onerror = () => reject(new Error("Cannot load " + src));
    image.src = src;
  });
}
function pixels(image) {
  // RGBA bytes of an image (exported PNGs carry no colour profile)
  const canvas = document.createElement("canvas");
  canvas.width = image.width; canvas.height = image.height;
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  ctx.drawImage(image, 0, 0);
  return ctx.getImageData(0, 0, image.width, image.height).data;
}
function twoChannels(rgba) {
  // Interleaved [R G] bytes: LUMINANCE_ALPHA texture data and CPU copy
  const out = new Uint8Array(rgba.length / 2);
  for (let i = 0, j = 0; i < rgba.length; i += 4, j += 2) { out[j] = rgba[i]; out[j + 1] = rgba[i + 1]; }
  return out;
}

// Sun: declination and the hour-angle offset of the current UTC time
// (low-precision solar ephemeris, ~0.1 deg)
function sunTerms(date) {
  const d = (date.getTime() - Date.UTC(2000, 0, 1, 12)) / 86400000;
  const g = (357.529 + 0.98560028 * d) * Math.PI / 180;
  const q = 280.459 + 0.98564736 * d;
  const L = (q + 1.915 * Math.sin(g) + 0.020 * Math.sin(2 * g)) * Math.PI / 180;
  const e = (23.439 - 0.00000036 * d) * Math.PI / 180;
  const decl = Math.asin(Math.sin(e) * Math.sin(L));
  const ra = Math.atan2(Math.cos(e) * Math.sin(L), Math.cos(L)) * 180 / Math.PI;
  const eot = ((q - ra) % 360 + 540) % 360 - 180;                       // deg
  const utcHours = date.getUTCHours() + date.getUTCMinutes() / 60 + date.getUTCSeconds() / 3600;
  return { decl, h0: (utcHours - 12) * 15 + eot };
}
function sunFade(x, y, sun) {
  // Night-only opacity at uv (x, y), as the shader: (limit + 1 - elevation) / 2, clamped
  const m = state.meta, r = Math.min(Math.max(Math.floor(y * m.height), 0), m.height - 1);
  const lon = m.lon_west + x * (m.lon_east - m.lon_west);
  const sinElevation = state.rowSin[r] * Math.sin(sun.decl) + state.rowCos[r] * Math.cos(sun.decl) * Math.cos((sun.h0 + lon) * Math.PI / 180);
  const elevation = Math.asin(Math.min(Math.max(sinElevation, -1), 1)) * 180 / Math.PI;
  return Math.min(Math.max((m.sun_limit_deg + 1 - elevation) / 2, 0), 1);
}

// ------------------------------------------------------------ time
function frameParts(frame) {
  const n = state.meta.times.length;
  const i = Math.max(0, Math.min(Math.floor(frame), n - 2));
  return { i, tau: Math.min(Math.max(frame - i, 0), 1) };
}
function frameDate(frame) {
  const { i, tau } = frameParts(frame);
  const t0 = state.meta.dates[i], t1 = state.meta.dates[i + 1] || t0;
  return new Date(t0.getTime() + tau * (t1.getTime() - t0.getTime()));
}
function intervalSeconds(i) {
  return state.meta.times.length > 1 ? (state.meta.dates[i + 1].getTime() - state.meta.dates[i].getTime()) / 1000 : 0;
}
function decodeColumn(code) {
  // Column density (birds km^-2) of a red code, also for a bilinear (non-integer) code as in the shader
  const m = state.meta;
  return code < 0.5 ? 0 : m.column_floor * (Math.exp((code - 1) / 254 * Math.log1p(m.column_max / m.column_floor)) - 1);
}
const currentWind = () => state.winds[Math.min(frameParts(state.frame).i, state.winds.length - 1)];

// --------------------------------------------------------- WebGL layer
const VERTEX = `
uniform mat4 u_matrix;
attribute vec2 a_pos;
varying vec2 v_m;
void main() { v_m = a_pos; gl_Position = u_matrix * vec4(a_pos, 0.0, 1.0); }`;

const FRAGMENT = `
precision highp float;
uniform sampler2D u_col0, u_col1, u_wind, u_lut;
uniform vec2 u_m0, u_m1;
uniform float u_tau, u_dt, u_morph, u_logmax, u_floor, u_climit, u_wind_limit, u_opacity;
uniform float u_night, u_sun_limit, u_decl, u_h0, u_lowfade;
varying vec2 v_m;
const float PI = 3.141592653589793;
const float EARTH = 40075016.686;
float column(float c) {
  float code = c * 255.0;
  return code < 0.5 ? 0.0 : u_floor * (exp((code - 1.0) / 254.0 * u_logmax) - 1.0);
}
vec2 windAt(vec2 uv) { vec4 t = texture2D(u_wind, uv); return vec2(t.r, t.a) * 2.0 * u_wind_limit - u_wind_limit; }
// Departure point of a parcel that moves s seconds with wind w (m/s, east/north)
vec2 depart(vec2 uv, vec2 w, float s, float k) { return uv - vec2(w.x, -w.y) * s * k / (u_m1 - u_m0); }
void main() {
  vec2 uv = (v_m - u_m0) / (u_m1 - u_m0);
  if (uv.x < 0.0 || uv.y < 0.0 || uv.x > 1.0 || uv.y > 1.0) discard;
  float mask = texture2D(u_col0, uv).a;
  if (mask < 0.02) discard;
  float y = (1.0 - 2.0 * v_m.y) * PI;
  float lat = atan(0.5 * (exp(y) - exp(-y)));
  float k = 1.0 / (EARTH * cos(lat));               // Mercator units per metre
  float value;
  if (u_morph > 0.5) {
    vec2 w = windAt(uv);
    float s0 = u_tau * u_dt;
    float s1 = -(1.0 - u_tau) * u_dt;
    vec2 uv0 = depart(uv, windAt(depart(uv, w, 0.5 * s0, k)), s0, k);
    vec2 uv1 = depart(uv, windAt(depart(uv, w, 0.5 * s1, k)), s1, k);
    value = (1.0 - u_tau) * column(texture2D(u_col0, uv0).r) + u_tau * column(texture2D(u_col1, uv1).r);
  } else {
    value = mix(column(texture2D(u_col0, uv).r), column(texture2D(u_col1, uv).r), u_tau);
  }
  float t = clamp(value / u_climit, 0.0, 1.0);
  vec3 color = texture2D(u_lut, vec2(t * 255.0 / 256.0 + 0.5 / 256.0, 0.5)).rgb;
  float alpha = u_opacity * smoothstep(0.0, 1.0, mask);
  if (u_lowfade > 0.5) alpha *= clamp(t * 8.0, 0.0, 1.0);
  if (u_night > 0.5) {
    float lon = v_m.x * 360.0 - 180.0;
    float elevation = degrees(asin(sin(lat) * sin(u_decl) + cos(lat) * cos(u_decl) * cos(radians(u_h0 + lon))));
    alpha *= clamp((u_sun_limit + 1.0 - elevation) / 2.0, 0.0, 1.0);
  }
  gl_FragColor = vec4(color * alpha, alpha);
}`;

function compile(gl, type, source) {
  const shader = gl.createShader(type);
  gl.shaderSource(shader, source);
  gl.compileShader(shader);
  if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(shader));
  return shader;
}
function texture(gl, format, width, height, data) {
  const tex = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, tex);
  gl.texImage2D(gl.TEXTURE_2D, 0, format, width, height, 0, format, gl.UNSIGNED_BYTE, data);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  return tex;
}

const birdLayer = {
  id: "birds", type: "custom", renderingMode: "2d",
  onAdd(_map, gl) {
    this.gl = gl;
    const program = gl.createProgram();
    gl.attachShader(program, compile(gl, gl.VERTEX_SHADER, VERTEX));
    gl.attachShader(program, compile(gl, gl.FRAGMENT_SHADER, FRAGMENT));
    gl.linkProgram(program);
    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(program));
    this.program = program;
    this.loc = {};
    for (const name of ["u_matrix", "u_col0", "u_col1", "u_wind", "u_lut", "u_m0", "u_m1", "u_tau", "u_dt", "u_morph",
      "u_logmax", "u_floor", "u_climit", "u_wind_limit", "u_opacity", "u_night", "u_sun_limit", "u_decl", "u_h0", "u_lowfade"]) {
      this.loc[name] = gl.getUniformLocation(program, name);
    }
    this.aPos = gl.getAttribLocation(program, "a_pos");
    // Own vertex array, so the map's bindings are not touched
    this.vaoExt = gl.createVertexArray ? null : gl.getExtension("OES_vertex_array_object");
    this.vao = gl.createVertexArray ? gl.createVertexArray() : this.vaoExt.createVertexArrayOES();
    this.bindVao(this.vao);
    this.buffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
    gl.enableVertexAttribArray(this.aPos);
    gl.vertexAttribPointer(this.aPos, 2, gl.FLOAT, false, 0, 0);
    this.bindVao(null);
    this.upload();
  },
  bindVao(vao) {
    if (this.gl.bindVertexArray) this.gl.bindVertexArray(vao); else this.vaoExt.bindVertexArrayOES(vao);
  },
  upload() {
    // Textures of the current night: columns and winds as LUMINANCE_ALPHA, colour map as RGB
    const gl = this.gl, m = state.meta;
    if (this.textures) for (const t of this.textures.all) gl.deleteTexture(t);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    const cols = state.cols.map((c) => texture(gl, gl.LUMINANCE_ALPHA, m.width, m.height, c));
    const winds = state.winds.map((w) => texture(gl, gl.LUMINANCE_ALPHA, m.width, m.height, w));
    const lut = texture(gl, gl.RGB, 256, 1, new Uint8Array(m.colormap.flat()));
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 4);
    this.textures = { cols, winds, lut, all: [...cols, ...winds, lut] };
    gl.bindBuffer(gl.ARRAY_BUFFER, this.buffer);
    const [x0, y0] = state.m0, x1 = x0 + state.span[0], y1 = y0 + state.span[1];
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([x0, y0, x1, y0, x0, y1, x1, y1]), gl.STATIC_DRAW);
  },
  render(gl, matrix) {
    if (!this.textures) return;
    const m = state.meta, L = this.loc, { i, tau } = frameParts(state.frame);
    const sun = sunTerms(frameDate(state.frame));
    gl.useProgram(this.program);
    this.bindVao(this.vao);
    gl.uniformMatrix4fv(L.u_matrix, false, matrix);
    const units = [[this.textures.cols[i], L.u_col0], [this.textures.cols[i + 1] || this.textures.cols[i], L.u_col1],
      [this.textures.winds[Math.min(i, this.textures.winds.length - 1)], L.u_wind], [this.textures.lut, L.u_lut]];
    units.forEach(([tex, loc], unit) => {
      gl.activeTexture(gl.TEXTURE0 + unit);
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.uniform1i(loc, unit);
    });
    gl.uniform2fv(L.u_m0, state.m0);
    gl.uniform2fv(L.u_m1, [state.m0[0] + state.span[0], state.m0[1] + state.span[1]]);
    gl.uniform1f(L.u_tau, tau);
    gl.uniform1f(L.u_dt, intervalSeconds(i));
    gl.uniform1f(L.u_morph, state.morph && this.textures.winds.length ? 1 : 0);
    gl.uniform1f(L.u_logmax, Math.log1p(m.column_max / m.column_floor));
    gl.uniform1f(L.u_floor, m.column_floor);
    gl.uniform1f(L.u_climit, state.climit);
    gl.uniform1f(L.u_wind_limit, m.wind_limit);
    gl.uniform1f(L.u_opacity, state.opacity);
    gl.uniform1f(L.u_night, state.night ? 1 : 0);
    gl.uniform1f(L.u_sun_limit, m.sun_limit_deg);
    gl.uniform1f(L.u_decl, sun.decl);
    gl.uniform1f(L.u_h0, sun.h0);
    gl.uniform1f(L.u_lowfade, state.lowfade ? 1 : 0);
    gl.enable(gl.BLEND);
    gl.blendFunc(gl.ONE, gl.ONE_MINUS_SRC_ALPHA);
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    this.bindVao(null);
  },
};

// ------------------------------------------------- CPU twin of the shader
function bilinear(bytes, x, y, ch) {
  // Channel ch (0 red, 1 green) of [R G] bytes at uv (x, y), filtered as the
  // GPU's LINEAR filter: pixel centres at (c + 0.5) / width, clamped at the edge
  const m = state.meta;
  const fx = Math.min(Math.max(x * m.width - 0.5, 0), m.width - 1.001);
  const fy = Math.min(Math.max(y * m.height - 0.5, 0), m.height - 1.001);
  const c = Math.floor(fx), r = Math.floor(fy), ax = fx - c, ay = fy - r;
  const p = 2 * (r * m.width + c) + ch, q = p + 2 * m.width;
  return (1 - ay) * ((1 - ax) * bytes[p] + ax * bytes[p + 2]) + ay * ((1 - ax) * bytes[q] + ax * bytes[q + 2]);
}
function windAt(wind, x, y, out) {
  // Birds' mean wind (m/s, east and north) at uv (x, y) from the [u v] bytes
  const m = state.meta, scale = 2 * m.wind_limit / 255;
  out[0] = bilinear(wind, x, y, 0) * scale - m.wind_limit;
  out[1] = bilinear(wind, x, y, 1) * scale - m.wind_limit;
  return out;
}
function displayedColumn(x, y) {
  // Column density (birds km^-2) that the shader draws at uv (x, y); NaN outside the data
  const { i, tau } = frameParts(state.frame);
  const col0 = state.cols[i], col1 = state.cols[i + 1] || col0;
  if (x < 0 || y < 0 || x > 1 || y > 1 || bilinear(col0, x, y, 1) < 0.02 * 255) return NaN;
  if (!(state.morph && state.winds.length)) {
    return (1 - tau) * decodeColumn(bilinear(col0, x, y, 0)) + tau * decodeColumn(bilinear(col1, x, y, 0));
  }
  const wind = currentWind(), dt = intervalSeconds(i), [sx, sy] = state.span;
  const k = mercPerMetre(state.m0[1] + y * sy);       // at the drawn point, as the shader
  const depart = (s) => {
    windAt(wind, x, y, W);
    windAt(wind, x - W[0] * 0.5 * s * k / sx, y + W[1] * 0.5 * s * k / sy, W);
    return [x - W[0] * s * k / sx, y + W[1] * s * k / sy];
  };
  const [x0, y0] = depart(tau * dt), [x1, y1] = depart(-(1 - tau) * dt);
  return (1 - tau) * decodeColumn(bilinear(col0, x0, y0, 0)) + tau * decodeColumn(bilinear(col1, x1, y1, 0));
}
function columnAt(x, y, frame) {
  // Column density of one exported frame at the nearest pixel; NaN outside the data
  const m = state.meta, col = state.cols[frame];
  const c = Math.min(Math.floor(x * m.width), m.width - 1), r = Math.min(Math.floor(y * m.height), m.height - 1);
  const p = r * m.width + c;
  return col[2 * p + 1] > 127 ? state.decodeTable[col[2 * p]] : NaN;
}

// ------------------------------------------------------ screen transform
// The image is regular in Web Mercator, so on an unpitched map (maxPitch 0)
// image uv maps to canvas pixels by an affine transform: origin at the
// north-west corner, columns along a, rows along b (from three projected corners)
const canvas = el("flow");
const ctx = canvas.getContext("2d");
function screenTransform() {
  const m = state.meta, ratio = window.devicePixelRatio || 1;
  const p = map.project([m.lon_west, m.lat_north]), pa = map.project([m.lon_east, m.lat_north]), pb = map.project([m.lon_west, m.lat_south]);
  const T = { ox: p.x * ratio, oy: p.y * ratio, ax: (pa.x - p.x) * ratio, ay: (pa.y - p.y) * ratio, bx: (pb.x - p.x) * ratio, by: (pb.y - p.y) * ratio };
  T.det = T.ax * T.by - T.bx * T.ay;
  return T;
}
function toUv(T, sx, sy, out) {
  const dx = sx - T.ox, dy = sy - T.oy;
  out[0] = (dx * T.by - dy * T.bx) / T.det;
  out[1] = (T.ax * dy - T.ay * dx) / T.det;
  return out;
}
function moveUv(x, y, u, v, s, out) {
  // Position after s seconds at wind (u, v) m/s, in image uv (y grows southward)
  const k = mercPerMetre(state.m0[1] + y * state.span[1]) * s;
  out[0] = x + u * k / state.span[0];
  out[1] = y - v * k / state.span[1];
  return out;
}

// ---------------------------------------------------------------- particles
// Particles move at the birds' own speed. When the clock advances by s
// simulated seconds (real time x playback speed, or a step of the slider or
// the arrow keys), each particle moves by wind x s, converted to map units
// as in the shader (two stages, wind at the half-way point), so particles
// and the advected density move together at any zoom and playback speed.
// Each particle is a streak back along its wind over TRAIL_SECONDS of
// flight, so the streak length also shows the speed; it fades in and out
// over FADE_SECONDS at the start and end of its life.
// Particles start where there are birds in the view (probability by column
// density, night only) and restart when they leave the birds, the night or
// the view. Their number follows the birds in the view: the Particles
// setting is the count at the night's maximum. When the nearest frame or
// the view changes, blocks of the view holding more particles than their
// share of the birds lose the excess, so the particles keep following the
// density (dusk spreading west, a zoom out) instead of staying where they
// started.
const TRAIL_SECONDS = 3600;
const FADE_SECONDS = 1200;
const LIFE_SECONDS = [2 * 3600, 5 * 3600];
const BLOCKS = 8;                                   // rebalance blocks per side of the view
const particles = { n: 0, active: 0, x: null, y: null, u: null, v: null, age: null, life: null, t: null,
  weights: null, cdf: null, cdfKey: null, totals: null, totalsKey: null, refill: 0, refillCdf: null, win: null, viewStale: true, warm: true };

function viewWindow(T) {
  // Image rows and columns under the canvas, clamped to the image
  const m = state.meta, uv = [0, 0], xs = [], ys = [];
  for (const [sx, sy] of [[0, 0], [canvas.width, 0], [0, canvas.height], [canvas.width, canvas.height]]) {
    toUv(T, sx, sy, uv); xs.push(uv[0]); ys.push(uv[1]);
  }
  const clamp = (v, n) => Math.min(Math.max(v, 0), n - 1);
  return {
    c0: clamp(Math.floor(Math.min(...xs) * m.width), m.width), c1: clamp(Math.ceil(Math.max(...xs) * m.width) - 1, m.width),
    r0: clamp(Math.floor(Math.min(...ys) * m.height), m.height), r1: clamp(Math.ceil(Math.max(...ys) * m.height) - 1, m.height),
    empty: Math.max(...xs) <= 0 || Math.min(...xs) >= 1 || Math.max(...ys) <= 0 || Math.min(...ys) >= 1,
  };
}
function birdWeights(frame, win, out) {
  // Column density of a frame over the view (zero where the sun is up, if
  // night only), row by row into out; returns the sum.
  // Separable sun test: sin(elevation) = a(row) + b(row) cos(h0 + lon(column))
  const m = state.meta, col = state.cols[frame], sun = sunTerms(m.dates[frame]);
  const limit = Math.sin((m.sun_limit_deg + 1) * Math.PI / 180);
  const cosH = state.colLon.map((lon) => Math.cos((sun.h0 + lon) * Math.PI / 180));
  let sum = 0;
  for (let r = win.r0, q = 0; r <= win.r1; r++) {
    const a = state.rowSin[r] * Math.sin(sun.decl), b = state.rowCos[r] * Math.cos(sun.decl);
    for (let c = win.c0; c <= win.c1; c++, q++) {
      const p = r * m.width + c;
      let w = col[2 * p + 1] > 127 ? state.decodeTable[col[2 * p]] : 0;
      if (w > 0 && state.night && a + b * cosH[c] > limit) w = 0;
      out[q] = w;
      sum += w;
    }
  }
  return sum;
}
function spawnTable(T) {
  // Birds of the nearest frame over the view: cumulative table to draw new
  // particles, totals of every frame to set the count, rebalance of the
  // blocks. The view is read again only after the map stops moving.
  const p = particles;
  if (p.viewStale || !p.win) { p.win = viewWindow(T); p.viewStale = false; }
  const m = state.meta, frame = Math.round(state.frame), win = p.win;
  const viewKey = [state.night, win.c0, win.c1, win.r0, win.r1, win.empty].join();
  const key = frame + "," + viewKey;
  if (p.cdfKey === key) return;
  const nc = win.c1 - win.c0 + 1, nr = win.r1 - win.r0 + 1, size = win.empty ? 0 : nc * nr;
  if (!p.weights || p.weights.length !== size) p.weights = new Float64Array(size);
  if (p.totalsKey !== viewKey) {
    // Birds in the view at every frame (the count follows them, relative to the night's maximum)
    p.totals = m.dates.map((_, f) => (size ? birdWeights(f, win, p.weights) : 0));
    p.totalsKey = viewKey;
  }
  const total = size ? birdWeights(frame, win, p.weights) : 0;
  const cdf = new Float64Array(size);
  for (let q = 0, sum = 0; q < size; q++) { sum += p.weights[q]; cdf[q] = sum; }
  p.cdf = cdf;
  p.cdfKey = key;
  const peak = Math.max(...p.totals);
  p.active = peak > 0 ? Math.round(p.n * total / peak) : 0;
  rebalance(win, total);
}
function rebalance(win, total) {
  // Blocks of the view with more particles than their share of the birds
  // (with a margin for chance) restart the excess; the restarted particles
  // reappear in the blocks short of particles (refill table, by density
  // within each block)
  const m = state.meta, p = particles;
  if (win.empty || !(total > 0) || !p.active) return;
  const nc = win.c1 - win.c0 + 1, nr = win.r1 - win.r0 + 1;
  const blockOf = (c, r) => Math.min(BLOCKS - 1, Math.floor((c - win.c0) / nc * BLOCKS)) * BLOCKS +
    Math.min(BLOCKS - 1, Math.floor((r - win.r0) / nr * BLOCKS));
  const share = new Float64Array(BLOCKS * BLOCKS), count = new Float64Array(BLOCKS * BLOCKS);
  for (let r = win.r0, q = 0; r <= win.r1; r++) for (let c = win.c0; c <= win.c1; c++, q++) share[blockOf(c, r)] += p.weights[q];
  const block = new Int32Array(p.active).fill(-1);
  for (let k = 0; k < p.active; k++) {
    if (Number.isNaN(p.x[k])) continue;
    const c = Math.floor(p.x[k] * m.width), r = Math.floor(p.y[k] * m.height);
    if (c < win.c0 || c > win.c1 || r < win.r0 || r > win.r1) continue;
    block[k] = blockOf(c, r);
    count[block[k]] += 1;
  }
  let restarted = 0;
  for (let k = 0; k < p.active; k++) {
    const b = block[k];
    if (b < 0) continue;
    const target = p.active * share[b] / total;
    if (count[b] > 1.2 * target + 2 && Math.random() > target / count[b]) { p.x[k] = NaN; restarted++; }
  }
  p.refill = 0;
  if (!restarted) return;
  const shortfall = share.map((w, b) => Math.max(0, p.active * w / total - count[b]) / Math.max(w, 1e-12));
  const refill = new Float64Array(p.weights.length);
  for (let r = win.r0, q = 0, sum = 0; r <= win.r1; r++) {
    for (let c = win.c0; c <= win.c1; c++, q++) { sum += p.weights[q] * shortfall[blockOf(c, r)]; refill[q] = sum; }
  }
  if (refill[refill.length - 1] > 0) { p.refillCdf = refill; p.refill = restarted; }
  p.warm = true;
}
function spawn(k, warm) {
  // New particle at a random bird of the view (or of the blocks short of
  // particles after a rebalance); warm particles start at a random age
  const m = state.meta, win = particles.win, refill = particles.refill > 0;
  const cdf = refill ? particles.refillCdf : particles.cdf;
  if (refill) particles.refill -= 1;
  const total = cdf.length ? cdf[cdf.length - 1] : 0;
  if (!(total > 0)) { particles.x[k] = NaN; return; }
  const target = Math.random() * total;
  let lo = 0, hi = cdf.length - 1;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (cdf[mid] < target) lo = mid + 1; else hi = mid; }
  const nc = win.c1 - win.c0 + 1, r = win.r0 + Math.floor(lo / nc), c = win.c0 + (lo % nc);
  particles.x[k] = (c + Math.random()) / m.width;
  particles.y[k] = (r + Math.random()) / m.height;
  particles.life[k] = LIFE_SECONDS[0] + Math.random() * (LIFE_SECONDS[1] - LIFE_SECONDS[0]);
  particles.age[k] = warm ? Math.random() * particles.life[k] : 0;
  windAt(currentWind(), particles.x[k], particles.y[k], W);
  particles.u[k] = W[0]; particles.v[k] = W[1];
}
function advanceParticles(s) {
  // Move every particle by s simulated seconds (midpoint rule, as the shader's two stages)
  const wind = currentWind(), p = particles, mid = [0, 0], out = [0, 0];
  for (let k = 0; k < p.active; k++) {
    if (Number.isNaN(p.x[k])) continue;
    windAt(wind, p.x[k], p.y[k], W);
    moveUv(p.x[k], p.y[k], W[0], W[1], 0.5 * s, mid);
    windAt(wind, mid[0], mid[1], W);
    // Mercator scale at the half-way point
    const k1 = mercPerMetre(state.m0[1] + mid[1] * state.span[1]) * s;
    out[0] = p.x[k] + W[0] * k1 / state.span[0];
    out[1] = p.y[k] - W[1] * k1 / state.span[1];
    p.x[k] = out[0]; p.y[k] = out[1]; p.u[k] = W[0]; p.v[k] = W[1];
    p.age[k] += Math.abs(s);
  }
}
function drawParticles(T) {
  const n = state.particleCount, ratio = window.devicePixelRatio || 1, p = particles;
  if (p.n !== n) {
    for (const key of ["x", "y", "u", "v", "age", "life"]) p[key] = new Float64Array(n);
    p.x.fill(NaN); p.n = n; p.warm = true; p.cdfKey = null;
  }
  spawnTable(T);
  // Advance by the simulated time since the last draw; a jump (slider, loop) restarts all
  const t = frameDate(state.frame).getTime() / 1000, s = p.t === null ? 0 : t - p.t;
  p.t = t;
  if (Math.abs(s) > 1.5 * intervalSeconds(frameParts(state.frame).i)) { p.x.fill(NaN); p.warm = true; }
  else if (s !== 0) advanceParticles(s);
  p.x.fill(NaN, p.active);                          // beyond the count of this frame
  const frame = Math.round(state.frame), sun = sunTerms(frameDate(state.frame));
  const margin = 20 * ratio, tail = [0, 0];
  // Paths by density level (3) and fade step (3), each with a faint tail half and a bright head half
  const paths = Array.from({ length: 18 }, () => new Path2D());
  for (let k = 0; k < p.active; k++) {
    if (Number.isNaN(p.x[k])) { spawn(k, p.warm); continue; }
    const x = p.x[k], y = p.y[k];
    const density = x < 0 || y < 0 || x >= 1 || y >= 1 ? NaN : columnAt(x, y, frame);
    if (!(density > 0) || p.age[k] > p.life[k] || (state.night && sunFade(x, y, sun) <= 0)) { spawn(k, false); continue; }
    const hx = T.ox + x * T.ax + y * T.bx, hy = T.oy + x * T.ay + y * T.by;
    if (hx < -margin || hy < -margin || hx > canvas.width + margin || hy > canvas.height + margin) { spawn(k, p.warm); continue; }
    moveUv(x, y, p.u[k], p.v[k], -TRAIL_SECONDS, tail);
    const tx = T.ox + tail[0] * T.ax + tail[1] * T.bx, ty = T.oy + tail[0] * T.ay + tail[1] * T.by;
    const mx = 0.5 * (hx + tx), my = 0.5 * (hy + ty);
    const level = density < state.climit * 0.1 ? 0 : density < state.climit * 0.4 ? 1 : 2;
    const fade = Math.min(p.age[k], p.life[k] - p.age[k]) / FADE_SECONDS;
    const step = fade < 1 / 3 ? 0 : fade < 2 / 3 ? 1 : 2;
    const b = 2 * (3 * level + step);
    paths[b].moveTo(tx, ty); paths[b].lineTo(mx, my);
    paths[b + 1].moveTo(mx, my); paths[b + 1].lineTo(hx, hy);
  }
  p.warm = false;
  ctx.lineWidth = 1.4 * ratio;
  ctx.lineCap = "round";
  [0.4, 0.7, 0.95].forEach((alpha, level) => [0.25, 0.6, 1].forEach((fade, step) => {
    const b = 2 * (3 * level + step);
    ctx.strokeStyle = `rgba(255,255,255,${0.35 * alpha * fade})`; ctx.stroke(paths[b]);
    ctx.strokeStyle = `rgba(255,255,255,${alpha * fade})`; ctx.stroke(paths[b + 1]);
  }));
}

// ---------------------------------------------------------------- arrows
// A regular screen grid of arrows. Each arrow points along the birds' mean
// wind of the current interval and its length is proportional to the speed
// (ARROW_PX_PER_MS CSS px per m/s, the same at every zoom; the legend shows
// ARROW_REF). Arrows are drawn where the displayed (morphed) column density
// is at least ARROW_MIN_DENSITY, more opaque where it is denser.
const ARROW_SPACING = 44;                 // CSS px between arrows
const ARROW_PX_PER_MS = 2.6;              // CSS px of arrow length per m/s (19 m/s, the fastest of us_east, ~ one spacing)
const ARROW_REF = 10;                     // m/s, legend arrow
const ARROW_MIN_DENSITY = 1;              // birds km^-2

function arrowPath(path, cx, cy, ex, ey, length, ratio) {
  // Shaft of the given length centred on (cx, cy) along unit (ex, ey), open head at the front
  const hx = cx + 0.5 * length * ex, hy = cy + 0.5 * length * ey;
  path.moveTo(cx - 0.5 * length * ex, cy - 0.5 * length * ey);
  path.lineTo(hx, hy);
  const h = Math.min(Math.max(0.35 * length, 3.5 * ratio), 9 * ratio), c = Math.cos(0.5), s = Math.sin(0.5);
  path.moveTo(hx - h * (ex * c - ey * s), hy - h * (ex * s + ey * c));
  path.lineTo(hx, hy);
  path.lineTo(hx - h * (ex * c + ey * s), hy - h * (-ex * s + ey * c));
}
function strokeArrows(context, paths, alphas, ratio) {
  // Dark halo under a white line, so the arrows read on every colour of the map
  context.lineCap = "round"; context.lineJoin = "round";
  paths.forEach((path, b) => {
    context.lineWidth = 3.4 * ratio; context.strokeStyle = `rgba(20,30,40,${0.55 * alphas[b]})`; context.stroke(path);
  });
  paths.forEach((path, b) => {
    context.lineWidth = 1.5 * ratio; context.strokeStyle = `rgba(255,255,255,${alphas[b]})`; context.stroke(path);
  });
}
function drawArrows(T) {
  const ratio = window.devicePixelRatio || 1, step = ARROW_SPACING * ratio, wind = currentWind();
  const sun = sunTerms(frameDate(state.frame)), uv = [0, 0], [spanX, spanY] = state.span;
  const alphas = [0.4, 0.6, 0.8, 1], paths = alphas.map(() => new Path2D());
  for (let sy = step / 2; sy < canvas.height; sy += step) {
    for (let sx = step / 2; sx < canvas.width; sx += step) {
      toUv(T, sx, sy, uv);
      const x = uv[0], y = uv[1], value = displayedColumn(x, y);
      if (!(value >= ARROW_MIN_DENSITY)) continue;
      let alpha = Math.min(1, 0.25 + 0.75 * value / (0.25 * state.climit));
      if (state.night) alpha *= sunFade(x, y, sun);
      if (alpha < 0.2) continue;
      windAt(wind, x, y, W);
      const speed = Math.hypot(W[0], W[1]);
      if (speed < 0.2) continue;
      // Screen direction of the wind (Mercator is conformal: a metre east and a metre north have the same screen length)
      const k = mercPerMetre(state.m0[1] + y * spanY), du = W[0] * k / spanX, dv = -W[1] * k / spanY;
      const ex = du * T.ax + dv * T.bx, ey = du * T.ay + dv * T.by, norm = Math.hypot(ex, ey);
      arrowPath(paths[Math.min(3, Math.round(alpha * 4) - 1)], sx, sy, ex / norm, ey / norm, speed * ARROW_PX_PER_MS * ratio, ratio);
    }
  }
  strokeArrows(ctx, paths, alphas, ratio);
}

// ---------------------------------------------------------------- flow overlay
// Redrawn after every map render (also while panning), so it stays on the map
function drawFlow() {
  if (!state.meta) return;
  const ratio = window.devicePixelRatio || 1;
  const width = Math.round(canvas.clientWidth * ratio), height = Math.round(canvas.clientHeight * ratio);
  if (canvas.width !== width || canvas.height !== height) { canvas.width = width; canvas.height = height; particles.viewStale = true; }
  ctx.clearRect(0, 0, width, height);
  if (state.flow === "particles") drawParticles(screenTransform());
  else if (state.flow === "arrows") drawArrows(screenTransform());
}
function drawFlowKey() {
  // Legend: a reference arrow or a particle streak, drawn on a low-density colour as on the map
  const key = el("flow-key"), c = key.getContext("2d"), ratio = window.devicePixelRatio || 1;
  key.width = Math.round(key.clientWidth * ratio); key.height = Math.round(key.clientHeight * ratio);
  const rgb = state.meta.colormap[40];
  c.fillStyle = `rgb(${rgb[0]},${rgb[1]},${rgb[2]})`;
  c.beginPath(); c.roundRect(0, 0, key.width, key.height, 3 * ratio); c.fill();
  const cx = key.width / 2, cy = key.height / 2;
  if (state.flow === "arrows") {
    const path = new Path2D();
    arrowPath(path, cx, cy, 1, 0, ARROW_REF * ARROW_PX_PER_MS * ratio, ratio);
    strokeArrows(c, [path], [1], ratio);
    el("flow-key-label").textContent = `${ARROW_REF} m/s (arrow length ∝ speed)`;
  } else if (state.flow === "particles") {
    const half = 0.38 * key.width;
    c.lineWidth = 1.4 * ratio; c.lineCap = "round";
    c.strokeStyle = "rgba(255,255,255,0.33)"; c.beginPath(); c.moveTo(cx - half, cy); c.lineTo(cx, cy); c.stroke();
    c.strokeStyle = "rgba(255,255,255,0.95)"; c.beginPath(); c.moveTo(cx, cy); c.lineTo(cx + half, cy); c.stroke();
    el("flow-key-label").textContent = "1 h of flight, moving at the birds' speed";
  }
}

// ---------------------------------------------------------------- runs
// Each run is an export folder (bmm4d.export_webmap); serve.py lists
// them with case, run, note and provenance, and ?data=<id> opens one
const escapeHtml = (text) => String(text).replace(/[&<>"]/g, (ch) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[ch]));
const shortDate = (iso) => new Date(iso).toLocaleString("en-GB", { day: "numeric", month: "short", year: "numeric", hour: "2-digit", minute: "2-digit" });
function fillRunPicker(selected) {
  // Runs grouped by root and case, the groups in the order of their newest run
  const groups = new Map();
  for (const r of state.runs) {
    const key = r.case ? `${r.root} · ${r.case}` : "";
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(r);
  }
  const option = (r) => `<option value="${escapeHtml(r.id)}" title="${escapeHtml(r.note || r.id)}"${r.id === selected ? " selected" : ""}>` +
    escapeHtml([r.run || r.id, r.created && shortDate(r.created)].filter(Boolean).join(" · ")) + "</option>";
  el("dataset").innerHTML = [...groups].map(([key, runs]) => (key ?
    `<optgroup label="${escapeHtml(key)}">${runs.map(option).join("")}</optgroup>` : runs.map(option).join(""))).join("");
}
function showRunInfo(id) {
  // Heading (case and run) and info lines (note, provenance, id) of a run from the server's index
  const r = state.runs.find((x) => x.id === id);
  if (config.title || !r || !r.case) return;
  el("heading").textContent = `${r.case} · ${r.run}`;
  el("heading").hidden = false;
  const made = [r.created && `exported ${shortDate(r.created)}`, r.job ? `cluster job ${r.job}` : r.host, r.commit && `commit ${r.commit}`];
  el("credit").textContent = [r.note, made.filter(Boolean).join(", "), r.id].filter(Boolean).join("\n");
  el("credit").hidden = false;
}

// ---------------------------------------------------------------- data
async function fetchNight(name) {
  // Meta and images of one exported night, decoded to [R G] bytes
  el("loading").hidden = false;
  const base = `data/${name}/`;
  const meta = await fetch(base + "meta.json").then((r) => r.json());
  meta.dates = meta.times.map((t) => new Date(t));
  const n = meta.times.length;
  const files = [...Array(n).keys()].map((i) => base + `col_${pad(i)}.png`)
    .concat([...Array(Math.max(n - 1, 0)).keys()].map((i) => base + `wind_${pad(i)}.png`));
  let done = 0;
  const images = await Promise.all(files.map((f) => loadImage(f).then((image) => {
    done += 1;
    el("loading-bar").style.width = `${100 * done / files.length}%`;
    el("loading-text").textContent = `Loading ${name}: ${done} / ${files.length}`;
    return image;
  })));
  return { meta, cols: images.slice(0, n).map((im) => twoChannels(pixels(im))),
    winds: images.slice(n).map((im) => twoChannels(pixels(im))) };
}
function showNight(night) {
  const meta = night.meta, n = meta.times.length;
  const switching = Boolean(state.meta);
  state.meta = meta;
  state.cols = night.cols;
  state.winds = night.winds;
  state.decodeTable = Float64Array.from({ length: 256 }, (_, code) => decodeColumn(code));
  // Image corner and span in Mercator units, latitude of each row, longitude of each column
  const nw = lib.MercatorCoordinate.fromLngLat([meta.lon_west, meta.lat_north]);
  const se = lib.MercatorCoordinate.fromLngLat([meta.lon_east, meta.lat_south]);
  state.m0 = [nw.x, nw.y]; state.span = [se.x - nw.x, se.y - nw.y];
  const rowLat = Float64Array.from({ length: meta.height }, (_, r) => latOfMercY(nw.y + (r + 0.5) / meta.height * (se.y - nw.y)));
  state.rowSin = rowLat.map((lat) => Math.sin(lat * Math.PI / 180));
  state.rowCos = rowLat.map((lat) => Math.cos(lat * Math.PI / 180));
  state.colLon = Float64Array.from({ length: meta.width }, (_, c) => meta.lon_west + (c + 0.5) / meta.width * (meta.lon_east - meta.lon_west));
  Object.assign(particles, { n: 0, t: null, cdfKey: null, totalsKey: null, viewStale: true });
  state.frame = Math.min(state.frame, n - 1);
  // Keep the colour limit when switching runs so they compare on one scale
  state.climit = switching ? Math.min(state.climit, meta.column_max) : meta.column_limit;
  el("time").max = String(n - 1);
  el("cmax").value = String(climitToSlider(state.climit));
  drawColorbar();
  drawFlowKey();
  updateTime();
  if (map.getLayer("birds")) birdLayer.upload(); else map.addLayer(birdLayer, layerAboveBirds());
  setRadars(meta.radars);
  if (!state.keepView) map.fitBounds([[meta.lon_west, meta.lat_south], [meta.lon_east, meta.lat_north]], { padding: 20, duration: 0 });
  state.keepView = true;
  el("loading").hidden = true;
  map.triggerRepaint();
}
function layerAboveBirds() {
  // The birds go under the first label layer (hosted styles) or under the
  // first layer marked "above-birds" (local style: borders and coasts)
  const layer = map.getStyle().layers.find((l) => l.type === "symbol" || (l.metadata && l.metadata["bmm:above-birds"]));
  return layer ? layer.id : undefined;
}
function setRadars(radars) {
  // A one-radar case has a single object in meta.json (MATLAB's jsonencode of a 1x1 struct)
  const data = { type: "FeatureCollection", features: [].concat(radars || []).map((r) => ({
    type: "Feature", properties: { name: r.name }, geometry: { type: "Point", coordinates: [r.lon, r.lat] } })) };
  if (map.getSource("radars")) { map.getSource("radars").setData(data); return; }
  map.addSource("radars", { type: "geojson", data });
  map.addLayer({ id: "radars", type: "circle", source: "radars", layout: { visibility: "none" },
    paint: { "circle-radius": 3, "circle-color": "#1d2730", "circle-stroke-color": "#fff", "circle-stroke-width": 1 } });
  // Names need the style's glyphs (the local basemap has none)
  if (!map.getStyle().glyphs) return;
  map.addLayer({ id: "radar-names", type: "symbol", source: "radars", minzoom: 5,
    layout: { visibility: "none", "text-field": ["get", "name"], "text-size": 11, "text-offset": [0, 0.9], "text-anchor": "top" },
    paint: { "text-color": "#1d2730", "text-halo-color": "#fff", "text-halo-width": 1 } });
}

// ---------------------------------------------------------------- UI
const climitToSlider = (v) => Math.log(v / 10) / Math.log(state.meta.column_max / 10);
const sliderToClimit = (s) => 10 * Math.pow(state.meta.column_max / 10, s);
const COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
function drawColorbar() {
  const bar = el("colorbar"), c = bar.getContext("2d"), image = c.createImageData(256, 1);
  state.meta.colormap.forEach((rgb, k) => image.data.set([rgb[0], rgb[1], rgb[2], 255], 4 * k));
  for (let row = 0; row < bar.height; row++) c.putImageData(image, 0, row);
  const fmt = (v) => (v >= 100 ? Math.round(v) : v.toFixed(v >= 10 ? 0 : 1)).toLocaleString();
  el("tick-1").textContent = fmt(state.climit / 2);
  el("tick-2").textContent = "≥ " + fmt(state.climit);
}
function updateTime() {
  const date = frameDate(state.frame);
  el("time").value = String(state.frame);
  el("time-label").textContent = date.toLocaleString("en-GB", {
    timeZone: "UTC", day: "numeric", month: "short", year: "numeric", hour: "2-digit", minute: "2-digit" }) + " UTC";
}
function setPlaying(playing) {
  state.playing = playing;
  el("play").textContent = playing ? "❚❚" : "▶";
  el("play").setAttribute("aria-label", playing ? "Pause" : "Play");
}
function setFlow(mode) {
  state.flow = ["particles", "arrows", "off"].includes(mode) ? mode : "particles";
  el("flow-mode").value = state.flow;
  el("particle-count-field").hidden = state.flow !== "particles";
  el("flow-legend").hidden = state.flow === "off";
  if (state.meta) drawFlowKey();
  if (map) map.triggerRepaint();
}
function bindUi() {
  el("play").addEventListener("click", () => setPlaying(!state.playing));
  el("time").addEventListener("input", (e) => { state.frame = Number(e.target.value); updateTime(); map.triggerRepaint(); });
  el("speed").addEventListener("change", (e) => { state.speed = Number(e.target.value); });
  el("cmax").addEventListener("input", (e) => { state.climit = sliderToClimit(Number(e.target.value)); drawColorbar(); map.triggerRepaint(); });
  el("opacity").addEventListener("input", (e) => { state.opacity = Number(e.target.value); map.triggerRepaint(); });
  el("morph").addEventListener("change", (e) => { state.morph = e.target.checked; map.triggerRepaint(); });
  el("night").addEventListener("change", (e) => { state.night = e.target.checked; particles.cdfKey = null; map.triggerRepaint(); });
  el("lowfade").addEventListener("change", (e) => { state.lowfade = e.target.checked; map.triggerRepaint(); });
  el("flow-mode").addEventListener("change", (e) => setFlow(e.target.value));
  el("particle-count").addEventListener("input", (e) => { state.particleCount = Number(e.target.value); map.triggerRepaint(); });
  el("radars").addEventListener("change", (e) => {
    for (const id of ["radars", "radar-names"]) if (map.getLayer(id)) map.setLayoutProperty(id, "visibility", e.target.checked ? "visible" : "none");
  });
  el("dataset").addEventListener("change", (e) => {
    const url = new URL(location.href); url.searchParams.set("data", e.target.value); history.replaceState(null, "", url);
    showRunInfo(e.target.value);
    fetchNight(e.target.value).then(showNight);
  });
  el("refresh").addEventListener("click", async () => {
    // Read the exports again (new runs, on the laptop or the cluster share)
    state.runs = await fetch("data/index.json?refresh=1").then((r) => r.json());
    fillRunPicker(el("dataset").value);
  });
  document.addEventListener("keydown", (e) => {
    if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") return;
    const n = state.meta ? state.meta.times.length : 1;
    if (e.code === "Space") { e.preventDefault(); setPlaying(!state.playing); }
    if (e.code === "ArrowRight" || e.code === "ArrowLeft") {
      const step = (e.code === "ArrowRight" ? 1 : -1) / 6;      // 5 min of a 30-min interval
      state.frame = Math.min(Math.max(state.frame + step, 0), n - 1); updateTime(); map.triggerRepaint();
    }
  });
  // The flow overlay follows every map frame; after a move, particles refill the new view
  map.on("render", drawFlow);
  map.on("moveend", () => { particles.viewStale = true; particles.warm = true; map.triggerRepaint(); });
  map.on("mousemove", (e) => {
    if (!state.meta) return;
    const p = lib.MercatorCoordinate.fromLngLat(e.lngLat);
    const x = (p.x - state.m0[0]) / state.span[0], y = (p.y - state.m0[1]) / state.span[1];
    const tip = el("tooltip"), value = displayedColumn(x, y);
    if (Number.isNaN(value)) { tip.hidden = true; return; }
    tip.hidden = false;
    tip.style.left = e.point.x + "px"; tip.style.top = e.point.y + "px";
    windAt(currentWind(), x, y, W);
    const heading = (Math.atan2(W[0], W[1]) * 180 / Math.PI + 360) % 360;
    tip.textContent = `${Math.round(value)} birds km⁻²` +
      (value >= ARROW_MIN_DENSITY ? ` · ${Math.hypot(W[0], W[1]).toFixed(0)} m/s to ${COMPASS[Math.round(heading / 22.5) % 16]}` : "");
  });
  map.on("mouseout", () => { el("tooltip").hidden = true; });
}

// ---------------------------------------------------------------- loop
function tick(now) {
  const dt = state.lastNow === null ? 0 : Math.min((now - state.lastNow) / 1000, 0.25);
  state.lastNow = now;
  if (state.meta && state.playing && state.meta.times.length > 1) {
    const n = state.meta.times.length, { i } = frameParts(state.frame);
    state.frame += dt * state.speed / intervalSeconds(i);
    if (state.frame > n - 1) state.frame = 0;
    updateTime();
    map.triggerRepaint();
  }
  requestAnimationFrame(tick);
}

// ---------------------------------------------------------------- start
async function start() {
  // Basemap: OpenFreeMap (default), Mapbox (?mapbox_token=), or the local
  // Natural Earth style shipped with the data (?basemap=local; no network)
  const token = params.get("mapbox_token");
  const basemap = config.basemap || params.get("basemap");
  let style = "https://tiles.openfreemap.org/styles/positron";
  const options = {};
  if (basemap === "local") {
    lib = window.maplibregl;
    style = await fetch("data/basemap.json").then((r) => r.json());
  } else if (token) {
    // Mapbox GL JS (same custom-layer API) with a Mapbox basemap, kept in Web Mercator
    const link = document.createElement("link");
    link.rel = "stylesheet"; link.href = "https://api.mapbox.com/mapbox-gl-js/v3.7.0/mapbox-gl.css";
    document.head.appendChild(link);
    await loadScript("https://api.mapbox.com/mapbox-gl-js/v3.7.0/mapbox-gl.js");
    lib = window.mapboxgl;
    lib.accessToken = token;
    style = "mapbox://styles/mapbox/light-v11";
    options.projection = "mercator";
  } else {
    lib = window.maplibregl;
  }
  // Runs: the server's index of exports (newest first), or the shared build's list
  state.runs = config.runs ? config.runs.map((id) => ({ id })) : await fetch("data/index.json").then((r) => r.json());
  const wanted = params.get("data");
  const run = state.runs.find((r) => r.id === wanted) || state.runs[0];
  if (!run) throw new Error("no web map export found; export an inversion with bmm4d.export_webmap (webmap/README.md)");
  fillRunPicker(run.id);
  el("dataset-field").hidden = Boolean(config.runs) && state.runs.length < 2;
  showRunInfo(run.id);
  if (config.title) { el("heading").textContent = config.title; el("heading").hidden = false; }
  if (config.credit) { el("credit").textContent = config.credit; el("credit").hidden = false; }
  const hash = config.hash !== false;
  state.keepView = hash && location.hash.length > 1;  // a shared view (#zoom/lat/lon) is kept
  map = new lib.Map({ container: "map", style, center: [-88, 38], zoom: 4, hash, maxPitch: 0,
    maxBounds: config.maxBounds, attributionControl: config.attribution === false ? false : { compact: true }, ...options });
  map.addControl(new lib.NavigationControl({ showCompass: false }), "top-right");
  map.addControl(new lib.ScaleControl({ unit: "metric" }), "bottom-right");
  setFlow(params.get("flow") || config.flow || state.flow);
  // The basemap and the night load together
  // (the birds go on the map once its style is ready, not after all basemap tiles)
  const styleReady = new Promise((resolve) => (map.isStyleLoaded() ? resolve() : map.once("style.load", resolve)));
  const [, night] = await Promise.all([styleReady, fetchNight(run.id)]);
  bindUi();
  showNight(night);
  requestAnimationFrame(tick);
}
start().catch((error) => {
  el("loading-text").textContent = "Error: " + error.message;
  console.error(error);
});
