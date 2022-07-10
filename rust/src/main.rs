

//Makes a CSR format sparse matrix out of a 3D laplacian
fn make_matrix(mx : usize, my : usize, mz : usize) -> (Vec<usize>,Vec<usize>,Vec<f64>){
    let mut cids  = Vec::<usize>::new();
    let mut offs  = Vec::<usize>::new();
    let mut vals  = Vec::<f64>::new();
    let id = |ix : usize, iy : usize, iz : usize| -> usize{
        iz + mz*(iy + my*ix)
    };

    let mut off : usize = 0;
    offs.push(off);

    for ix in 0..mx{
        for iy in 0..my{
            for iz in 0..mz{
                {
                    let mut add_nonzero = |ix : usize, iy : usize, iz : usize, val : f64|{
                        vals.push(val);
                        cids.push(id(ix,iy,iz));
                        off+=1;
                    };

                    add_nonzero(ix,iy,iz,6.0);
                    if ix>0{
                        add_nonzero(ix-1,iy,iz,-1.0);
                    }
                    if ix<mx-1{
                        add_nonzero(ix+1,iy,iz,-1.0);
                    }
                    if iy>0{
                        add_nonzero(ix,iy-1,iz,-1.0);
                    }
                    if iy<my-1{
                        add_nonzero(ix,iy+1,iz,-1.0);
                    }
                    if iz>0{
                        add_nonzero(ix,iy,iz-1,-1.0);
                    }
                    if iz<mz-1{
                        add_nonzero(ix,iy,iz+1,-1.0);
                    }
                }
                offs.push(off);
            }
        }
    }

    (cids,offs,vals)
}


fn csr_eval(cids : &[usize], offs : &[usize], vals : &[f64], x : &[f64], y : &mut [f64]){
    for (yi,(&beg,&end)) in y.iter_mut().zip(offs.iter().zip(offs.iter().skip(1))){
        *yi = 0.0;
        for (val,&ri) in vals[beg..end].iter().zip(cids[beg..end].iter()){
            *yi += val * x[ri];
        }
    }
}

fn dot(x : &[f64], y : &[f64]) -> f64{
    x.iter().zip(y.iter()).map(|(x,y)| x*y).sum()
}


fn cgsolve(cids : &[usize], offs : &[usize], vals : &[f64], b : &[f64], r : &mut [f64], p : &mut [f64], ap : &mut[f64], x : &mut [f64]) -> f64{
    let tol = 1e-8;
    //Set r = b - A*x 
    csr_eval(cids,offs,vals,x,r);
    for (ri,bi) in r.iter_mut().zip(b.iter()){
        *ri = *bi - *ri;
    }
    //Set p = r
    for (pi,ri) in p.iter_mut().zip(r.iter()){
        *pi = *ri;
    }
    //Set rho = r^T*r
    let mut rho = dot(r,r);
    while rho.sqrt() > tol{
        csr_eval(cids,offs,vals,p,ap);
        let rdotr = rho;
        let pdotap = dot(p,ap);
        let alpha = rdotr/pdotap;
        for (xi,&pi) in x.iter_mut().zip(p.iter()){
            *xi+=alpha*pi;
        }
        for (ri,&api) in r.iter_mut().zip(ap.iter()){
            *ri-=alpha*api;
        }
        let rpdotrp = dot(r,r);
        let beta = rpdotrp/rdotr;
        for (pi, &ri) in p.iter_mut().zip(r.iter()){
            *pi=ri+beta*(*pi);
        }
        rho = dot(r,r);
    }
    rho.sqrt()
}



fn main() {
    use std::env;
    let args : Vec<String> = env::args().collect();
    if args.len()==4{
        let mx = args[1].parse::<usize>().unwrap();
        let my = args[2].parse::<usize>().unwrap();
        let mz = args[3].parse::<usize>().unwrap();
        let m = mx*my*mz;
        let (cids,offs,vals) = make_matrix(mx,my,mz);
        let mut x = vec![0.0;m];
        let mut b = vec![1.0;m];
        let mut r = vec![0.0;m];
        let mut p = vec![0.0;m];
        let mut ap = vec![0.0;m];
        let res = cgsolve(&cids,&offs,&vals,&b,&mut r,&mut p,&mut ap,&mut x);
        println!("res : {}",res);
    }
    else{
        println!("Usage: ./main mx my mz");
    }
}
