use plotly::{Plot, Scatter};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;

// Model Params
const N: usize = 16;
const SEED: u64 = 42;
const J: f32 = 1.0; // Coupling constant
const TEMP_MAX: f32 = 5.0;
const DT: f32 = 0.05; // Temperature Difference
const N_SWEEP: usize = 2000; // MCMC Sweep
const N_BINS: usize = 1;
const NUM_T: usize = ((TEMP_MAX - 0.1) / DT + 1.0) as usize;

// Initialize a N X N matrix
fn init_lattice() -> [[f32; N]; N] {
    // let mut rng = rand::thread_rng();
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut matrix: [[f32; N]; N] = [[0.0; N]; N];

    for i in 0..N {
        for j in 0..N {
            let rng: i64 = rng.gen();
            if (rng % 2) == 0 {
                matrix[i][j] = 1.0;
            } else {
                matrix[i][j] = -1.0;
            }
        }
    }
    matrix
}

fn get_delta_e(lattice: &mut [[f32; N]; N], i: usize, j: usize) -> f32 {
    let curr_site: f32 = lattice[i][j];
    let adj_sites: f32 = lattice[(i + 1) % (N)][j]
        + lattice[(i + N - 1) % (N)][j]
        + lattice[i][(j + 1) % (N)]
        + lattice[i][(j + N - 1) % (N)];
    let delta_e: f32 = 2.0 * J * curr_site * adj_sites;
    delta_e
}

fn flip_lattice(lattice: &mut [[f32; N]; N], curr_temp: f32) {
    let mut rng = rand::thread_rng();
    for _i in 0..N.pow(2) {
        let nx: usize = rng.gen_range(0..=(N - 1));
        let ny: usize = rng.gen_range(0..=(N - 1));
        let delta_e: f32 = get_delta_e(lattice, nx, ny);
        let rand_val: f32 = rng.gen();
        if (-1.0 * delta_e * curr_temp.powf(-1.0)).exp() > rand_val {
            lattice[nx][ny] = -1.0 * lattice[nx][ny];
        }
    }
}

fn cal_magnetization(
    lattice: &mut [[f32; N]; N],
    magnetization_mat: &mut [[f32; N_BINS * N_SWEEP]; NUM_T],
    temp_i: usize,
    s: usize,
) {
    let mut m: f32 = 0.0;
    for i in 0..N {
        for j in 0..N {
            m = m + lattice[i][j];
        }
    }
    magnetization_mat[temp_i][s] = m / (N.pow(2) as f32)
}

fn metropolis_sampling(
    lattice: &mut [[f32; N]; N],
    magnetization_mat: &mut [[f32; N_BINS * N_SWEEP]; NUM_T],
    temp_array: [f32; NUM_T],
    temp_i: usize,
) {
    let curr_temp: f32 = temp_array[temp_i];
    for s in 0..N_SWEEP {
        flip_lattice(lattice, curr_temp);
        cal_magnetization(lattice, magnetization_mat, temp_i, s as usize);
    }
    println!(
        "Complete Sampling in current temperature {:.2}",
        temp_array[temp_i]
    );
}

// Handy function for us to check matrix
fn print_matrix(matrix: Vec<Vec<f32>>) {
    for row in &matrix {
        print!("[");
        for &element in row {
            print!("{:.1} ", element);
        }
        println!("],");
    }
}

// Test random number generatto
fn main() {
    let mut lattice: [[f32; N]; N] = init_lattice();
    let mut t_array: [f32; NUM_T] = [0.0; NUM_T];
    let mut magnetization_mat: [[f32; N_BINS * N_SWEEP]; NUM_T] = [[0.0; N_BINS * N_SWEEP]; NUM_T];
    let mut aver_m: [f32; NUM_T] = [0.0; NUM_T];

    // Performe MH Sampling
    for i in 0..NUM_T {
        t_array[i] = 0.1 + (i as f32) * DT;
        metropolis_sampling(&mut lattice, &mut magnetization_mat, t_array, i);
    }

    // Compute AVG magnetization
    for k in 0..NUM_T {
        let mut m_avg: f32 = 0.0;
        for s in 0..N_SWEEP * N_BINS {
            m_avg = m_avg + magnetization_mat[k][s];
        }
        aver_m[k] = m_avg / (N_SWEEP as f32);
    }

    // Print result
    for i in 0..NUM_T {
        println!(
            "temp {:.2}, abs_magetization{:.2}",
            t_array[i],
            aver_m[i].abs()
        )
    }

    let mut abs_aver_m: Vec<f32> = vec![0.0; NUM_T];
    for i in 0..NUM_T {
        abs_aver_m[i] = aver_m[i].abs()
    }

    let mut plot = Plot::new();
    let trace = Scatter::new(t_array.to_vec(), abs_aver_m.to_vec());
    plot.add_trace(trace);
    plot.show()
}

// fn main() {
//     let mut lattice: Vec<Vec<f32>> = init_lattice();
//     let curr_temp: f32 = 0.05;
//     print_matrix(lattice.clone());
//     // perform flip
//     flip_lattice(&mut lattice, curr_temp);
//     println!("Flip matric");
//     print_matrix(lattice.clone());
//
// }
