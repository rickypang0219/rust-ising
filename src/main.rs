use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use plotly::{Plot, Scatter};


// Model Params
const N: usize = 32;
const SEED: u64 = 42;
const J:f32 = 1.0; // Coupling constant
const TEMP_MAX: f32 = 5.0;
const DT: f32 = 0.01; // Temperature Difference
const N_SWEEP: usize =  1000; // MCMC Sweep
const N_BINS:usize = 5;
const NUM_T:usize = ( (TEMP_MAX - 0.1) / DT + 1.0) as usize;



// Initialize a N X N matrix
fn init_lattice() -> Vec<Vec<f32>> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let mut matrix: Vec<Vec<f32>> = vec![vec![0.0;N];N];

    for i in 0..N {
        for j in 0..N {
            let rng:u64 = rng.gen();
            if (rng % 2) == 1 {
                matrix[i][j] = 1.0;
            } else{
                matrix[i][j] = -1.0;
            }
        }
    }
    return matrix;
}



fn get_delta_e(
    lattice:&mut Vec<Vec<f32>>,
    i:usize,
    j:usize) -> f32 {
    let curr_site: f32 = lattice[i][j];
    let adj_sites: f32= lattice[(i + 1)%(N)][j] + lattice[(i + N - 1)%(N)][j] + lattice[i][(j+1)%(N)] + lattice[i][(j + N - 1)%(N)];
    let delta_e: f32 =  2.0 * J * curr_site * adj_sites ;
    return delta_e;
}


fn flip_lattice(
    lattice: &mut Vec<Vec<f32>>,
    curr_temp:f32){
    let mut rng = rand::thread_rng();
    for _i in 0..N.pow(2) {
        let nx: usize = rng.gen_range(0..=(N - 1)) ;
        let ny: usize = rng.gen_range(0..=(N -1 )) ;
        let delta_e : f32 = get_delta_e(
            lattice,
            nx, ny);
        let rand_val: f32 = rng.gen::<f32>();
        if (-1.0 *  delta_e / curr_temp).exp() > rand_val {
            lattice[nx][ny] = -1.0 * lattice[nx][ny];
        }
    }
}


fn cal_magnetization(
    lattice:&mut Vec<Vec<f32>>,
    magnetization_mat: &mut Vec<Vec<f32>>,
    temp_i:usize,
    s:usize) {
    let mut  m:f32 = 0.0;
    for i in 0..N {
        for j in 0..N {
            m = m + lattice[i][j];
        }
    }
    magnetization_mat[temp_i][s] = m / (N.pow(2) as f32)
}




fn metropolis_sampling(
    lattice: &mut Vec<Vec<f32>>,
    magnetization_mat:&mut Vec<Vec<f32>>,
    temp_array:Vec<f32>,
    temp_i:usize) {
    let curr_temp:f32 = temp_array[temp_i];
    for s in 0..N_SWEEP {
        flip_lattice(
            lattice,
            curr_temp);
        cal_magnetization(
            lattice,
            magnetization_mat,
            temp_i,
            s as usize);
    }
    println!("Complete Sampling in current temperature {:.2}", temp_array[temp_i]);

}



// Handy function for us to check matrix
fn print_matrix(matrix : Vec<Vec<f32>> ) {
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
    let mut lattice: Vec<Vec<f32>> = init_lattice();
    let mut t_array : Vec<f32> = vec![0.0;NUM_T];
    let mut magnetization_mat: Vec<Vec<f32>> = vec![vec![0.0;N_BINS * N_SWEEP];NUM_T];
    let mut aver_m:Vec<f32> = vec![0.0;NUM_T];

    // Performe MH Sampling
    for i in 0..NUM_T {
        t_array[i] = 0.1 + (i as f32)  * DT ;
        metropolis_sampling(
            &mut lattice,
            &mut magnetization_mat,
            t_array.clone(),
            i);
    }

    // Compute AVG magnetization
    for k in 0..NUM_T {
        let mut m_avg: f32 = 0.0;
        for s in 0..N_SWEEP {
            m_avg = m_avg + magnetization_mat[k][s];
        }
        aver_m[k] = m_avg / (N_SWEEP as f32);
    }

    // Print result
    for i in 0..NUM_T {
        println!("temp {:.2}, abs_magetization{:.2}", t_array[i], aver_m[i].abs() )
    }

    let mut abs_aver_m: Vec<f32> = vec![0.0;NUM_T];
    for i in 0..NUM_T {
        abs_aver_m[i] = aver_m[i].abs()
    }

    let mut plot = Plot::new();
    let trace = Scatter::new(t_array.clone(), abs_aver_m.clone());
    plot.add_trace(trace);
    plot.show()

}



// fn main() {
//     let mut lattice: Vec<Vec<f64>> = init_lattice();
//     let curr_temp: f64 = 0.05;
//     print_matrix(lattice.clone());
//     // perform flip
//     flip_lattice(&mut lattice, curr_temp);
//     println!("Flip matric");
//     print_matrix(lattice.clone());
//
// }
