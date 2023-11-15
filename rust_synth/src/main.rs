use std::fs::File;
use std::io::Read;

use hound;

fn read_bytes_f64_file(filename: &str) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let mut file = File::open(filename)?;
    let mut buf: Vec<u8> = Vec::new();
    let _ = file.read_to_end(&mut buf)?;
    let float_buf: Vec<f64> = buf
        .chunks(8)
        .map(|chunk| {
            let mut bytes_array: [u8; 8] = [0; 8];
            bytes_array.copy_from_slice(chunk);
            f64::from_le_bytes(bytes_array)
        })
        .collect();
    Ok(float_buf)
}

fn ode(
    t: usize,
    v: [f64; 6],
    pi: &mut Vec<f64>,
    pb: &mut Vec<f64>,
    alpha: &Vec<f64>,
    beta: &Vec<f64>,
    gamma: f64,
    envelope: &Vec<f64>,
    c: f64,
    L: f64,
    r: f64,
    Ch: f64,
    MG: f64,
    MB: f64,
    RB: f64,
    Rh: f64,
    dt: f64,
    ovfs: usize,
) -> [f64; 6] {
    let [x, y, pout, i1, i2, i3] = v;
    let dv0 = y;
    let dv1 = (-alpha[t / ovfs] - beta[t / ovfs] * x - x * x * x + x * x) * gamma * gamma
        - (x * x * y + x * y) * gamma;

    let pbold = pb[t];
    let delayed_t = if t > ((L / c / dt) as usize) {t - ((L / c / dt) as usize)} else {0};
    pi[t] = (0.5 * envelope[t / ovfs]) * dv1 + pb[delayed_t];
    pb[t] = -r * pi[delayed_t];
    let pout = (1.0 - r) * pi[delayed_t];

    let dv2 = (pb[t] - pbold) / dt;
    let dv3 = i2;
    let dv4 = -(1.0 / Ch / MG) * i1 - Rh * (1.0 / MB + 1.0 / MG) * i2
        + (1.0 / MG / Ch + Rh * RB / MG / MB) * i3
        + (1.0 / MG) * dv2
        + (Rh * RB / MG / MB) * pout;
    let dv5 = -(MG / MB) * i2 - (Rh / MB) * i3 + (1.0 / MB) * pout;
    [dv0, dv1, dv2, dv3, dv4, dv5]
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let alpha = read_bytes_f64_file("data/alpha.raw")?;
    let beta = read_bytes_f64_file("data/beta.raw")?;
    let envelope = read_bytes_f64_file("data/envelope.raw")?;

    let fs = 44100;
    let gamma = 40000.0;
    let ovfs = 20;

    let mut t = 0;
    let tmax = alpha.len() * ovfs - 1;
    let dt = 1.0 / (ovfs * fs) as f64;

    let mut pi: Vec<f64> = vec![0.0; tmax];
    let mut pb: Vec<f64> = vec![0.0; tmax];
    let mut out: Vec<f64> = vec![0.0; alpha.len()];

    let mut v = [1e-4 * 1e2, 1e-4 * 1e1, 1e-4, 1e-4, 1e-4, 1e-4];
    let (c, L, r, Ch, MG, MB, RB, Rh) = (
        343.0, 0.025, 0.65, 1.43e-10, 20.0, 10000.0, 5000000.0, 24000.0,
    );

    for t in 0..tmax {
        let k1 = ode(
            t, v, &mut pi, &mut pb, &alpha, &beta, gamma, &envelope, c, L, r, Ch, MG, MB, RB, Rh,
            dt, ovfs,
        );

        let v2_vec: Vec<f64> = v
            .iter()
            .zip(k1.iter())
            .map(|(v_, k1_)| v_ + dt / 2.0 * k1_)
            .collect();
        let mut v2: [f64; 6] = [0.0; 6];
        v2.copy_from_slice(&v2_vec);
        let k2 = ode(
            t, v2, &mut pi, &mut pb, &alpha, &beta, gamma, &envelope, c, L, r, Ch, MG, MB, RB, Rh,
            dt, ovfs,
        );

        let v3_vec: Vec<f64> = v
            .iter()
            .zip(k2.iter())
            .map(|(v_, k2_)| v_ + dt / 2.0 * k2_)
            .collect();
        let mut v3: [f64; 6] = [0.0; 6];
        v3.copy_from_slice(&v3_vec);
        let k3 = ode(
            t, v3, &mut pi, &mut pb, &alpha, &beta, gamma, &envelope, c, L, r, Ch, MG, MB, RB, Rh,
            dt, ovfs,
        );

        let v4_vec: Vec<f64> = v
            .iter()
            .zip(k3.iter())
            .map(|(v_, k3_)| v_ + dt * k3_)
            .collect();
        let mut v4: [f64; 6] = [0.0; 6];
        v4.copy_from_slice(&v4_vec);
        let k4 = ode(
            t, v4, &mut pi, &mut pb, &alpha, &beta, gamma, &envelope, c, L, r, Ch, MG, MB, RB, Rh,
            dt, ovfs,
        );

        v = [
            v[0] + dt * (2.0 * (k2[0] + k3[0]) + k1[0] + k4[0]) / 6.0,
            v[1] + dt * (2.0 * (k2[1] + k3[1]) + k1[1] + k4[1]) / 6.0,
            v[2] + dt * (2.0 * (k2[2] + k3[2]) + k1[2] + k4[2]) / 6.0,
            v[3] + dt * (2.0 * (k2[3] + k3[3]) + k1[3] + k4[3]) / 6.0,
            v[4] + dt * (2.0 * (k2[4] + k3[4]) + k1[4] + k4[4]) / 6.0,
            v[5] + dt * (2.0 * (k2[5] + k3[5]) + k1[5] + k4[5]) / 6.0,
        ];

        out[t / ovfs] = RB * v[5];
    }

    let max_abs_out: f64 = out.iter().map(|x| x.abs()).fold(f64::NEG_INFINITY, |a, b| a.max(b));
    let mut i16_out: Vec<i16> = vec![0; out.len()];
    for i in 0..out.len() {
        i16_out[i] = (out[i] / max_abs_out * i16::MAX as f64) as i16;
    }

    let spec = hound::WavSpec {
        channels: 1,           // モノラル
        sample_rate: 44100,    // サンプルレート
        bits_per_sample: 16,   // サンプルあたりのビット数
        sample_format: hound::SampleFormat::Int, // サンプルフォーマット
    };
    let mut writer = hound::WavWriter::create("out/out.wav", spec)?;
    for sample in i16_out.iter() {
        writer.write_sample(*sample)?;
    }

    Ok(())
}
