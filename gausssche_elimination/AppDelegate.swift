//
//  AppDelegate.swift
//  gausssche_elimination
//
//  Created by Peter Weilnböck on 24/07/16.
//  Copyright © 2016 peter.weilnboeck. All rights reserved.
//

import Cocoa

@NSApplicationMain
class AppDelegate: NSObject, NSApplicationDelegate {

    @IBOutlet weak var window: NSWindow!


    func applicationDidFinishLaunching(_ aNotification: Notification) {
        // Insert code here to initialize your application
        print("Hello, world!")
        let L = 16
        var M = Array<Array<Double>>()
        var p = Array(repeating: Int(), count: L)
        
        //Gauss'sche Elimination mit Pivoting:
        for _ in 0...(L-1){
            M.append(Array(repeating: Double(), count: L))
        }
        let us2 = 1.0/sqrt(2.0)
        
        var m01,m02,m03,m04,m05,m06,m07,m08,m09,m10,m11,m12,m13,m14,m15,m16,vB: Array<Double>
        var x,B: Array<Array<Double>>
        
        m01 = [0,1,1,0,0,0,us2,0,0,0,0,0,0,0,0,0]
        m02 = [1,0,0,0,0,0,-us2,0,0,0,0,0,0,0,0,0]
        m03 = [0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0]
        m04 = [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0]
        m05 = [0,0,0,1,1,0,0,0,-us2,0,us2,0,0,0,0,0]
        m06 = [0,0,0,0,0,0,0,0,-us2,-1,-us2,0,0,0,0,0]
        m07 = [0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0]
        m08 = [0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0]
        m09 = [0,0,0,0,0,1,0,0,0,0,0,0,-us2,0,0,0]
        m10 = [0,0,0,0,0,0,0,0,0,0,0,0,-us2,0,0,1]
        m11 = [0,0,0,0,0,0,us2,0,-us2,0,0,0,0,1,0,0]
        m12 = [0,0,0,0,0,0,-us2,-1,-us2,0,0,0,0,0,0,0]
        m13 = [0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0]
        m14 = [0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0]
        m15 = [0,0,0,0,0,0,0,0,0,0,us2,0,-us2,0,1,0]
        m16 = [0,0,0,0,0,0,0,0,0,0,-us2,-1,-us2,0,0,0]
        M = [m01,m02,m03,m04,m05,m06,m07,m08,m09,m10,m11,m12,m13,m14,m15,m16]
        
        
        let f = 1.0
        let f1 = f/2.0
        let f2 = f/6.0
        let f3 = f/2.0
        vB = [0,0,0,f1,0,f2,0,f3,0,0,0,0,0,0,0,0]
        
        for row in 0...(L-1){
            for column in 0...(L-1){
                print("Row \(row) column \(column) is \(M[row][column])")
            }
        }
        
        (M, p) = do_gauss_elemimination(M, L: L)
        
        B = [vB]
        B = transpose(B)
        x = resubstitute(B, M: M, p: p, L: L, LB: 1)
        dump(x)
        
        //Rücksubstitution bei Gauss'scher Eliminierung mit Pivoting
        //Initialisierung:
//        let LB = 3
//        var B = Array<Array<Double>>()
//        B = [[1,0,0],[0,1,0],[0,0,1]]
//        
//        let X = resubstitute(B, M: M, p: p, L: L, LB: LB)
//        dump(X)
        
        let invM = calculate_inverse_Matrix(M, p: p, L: L)
//        dump(invM)
        
        let detM = calculate_Determinante(M, p: p, L: L)
        
        print("Determinante von M =\(detM)")
    }
    
    func do_gauss_elemimination(_ inputM: Array<Array<Double>>, L: Int) -> (returnM: Array<Array<Double>>, returnP: Array<Int>){
        var p = Array(repeating: Int(), count: L)
        var M = inputM
        for n in 0...(L-1){
            p[n] = n
        }
        
        var piv = 0.0
        var ipiv = 0
        var merk = 0
        //Masterloop
        for n in 0...(L-2){
//            dump(M)
//            dump(p)
            //Suche das Pivot Element
            piv = 0.0
            for i in n...(L-1){
                if abs(M[p[i]][n])>piv {
                    piv = abs(M[p[i]][n])
                    ipiv = i
                }
            }
            //Virtuelle Vertauschung
            merk = p[n]
            p[n] = p[ipiv]
            p[ipiv] = merk
            //Gauss'sche Elimination
            var f = 0.0
            for i in (n+1)...(L-1){
                f = (M[p[i]][n])/(M[p[n]][n])
                //                M[p[i]][n] = f // Änderung hier - fehler in den Übungsaufgaben?
                for j in (n)...(L-1){
                    M[p[i]][j] = M[p[i]][j] - f * M[p[n]][j]
                }
            }
            
            
        }
        //        for row in 0...(L-1){
        //            for column in 0...(L-1){
        //                print("Row \(row) column \(column) is \(M[row][column])")
        //            }
        //        }
//        dump(M)
//        dump(p)
        return(M, p)
    }
    
    func resubstitute(_ inputB: Array<Array<Double>>, M: Array<Array<Double>>, p: Array<Int>, L: Int, LB: Int) -> Array<Array<Double>>{
        
        var B = inputB
        //Eliminationsschritte für die Rechte Seite B:
        for n in 0...(L-2){
            for i in (n+1)...(L-1){
                for j in 0...(LB-1){
                    B[p[i]][j] = B[p[i]][j] - M[p[i]][n] * B[p[n]][j]
                }
            }
        }
        //masterloop über Spalten von B:
        var x = Array<Array<Double>>()
//        x = [[1,0,0],[0,1,0],[0,0,1]]
        x = B
        for ib in 0...(LB-1){
            for n in stride(from: (L-1), to: 0, by: -1){
                x[n][ib] = B[p[n]][ib]
                if n<(L-1) {
                    for j in (n+1)...(L-1){
                        x[n][ib] = x[n][ib] - M[p[n]][j] * x[j][ib]
                    }
                }else{
                    x[n][ib] = x[n][ib] - M[p[n]][L-1] * x[L-1][ib]
                }
                x[n][ib] = x[n][ib]/M[p[n]][n]
            }
        }
        return x
    }
    
    func calculate_inverse_Matrix(_ inputM: Array<Array<Double>>, p: Array<Int>, L: Int) -> Array<Array<Double>>{
        var B = Array<Array<Double>>()
        for _ in 0...(L-1){
            B.append(Array(repeating: Double(), count: L))
        }
        for row in 0...(L-1){
            for column in 0...(L-1){
                if row == column {
                    B[row][column] = 1.0
                }else{
                    B[row][column] = 0.0
                }
            }
        }
        let invM = resubstitute(B, M: inputM, p: p, L: L, LB: L)
        
        return invM
    }
    
    func calculate_Determinante(_ inputM: Array<Array<Double>>, p: Array<Int>, L: Int) -> Double{
        var M = inputM
        var pf = p
        //Vorzeichen für Determinante
        var detM = 1.0
        for n in 0...(L-1){
            detM = detM * M[pf[n]][n]
//            print("Determinante von M =\(detM), n=\(n), M=\(M[p[n]][n])")
        }
        var s = 1.0
        var merk = 0
        for i in 0...(L-1){
            while i != pf[i] {
                merk = pf[i]
                pf[i] = pf[merk]
                pf[merk] = merk
                s = -s
            }
        }
        detM = s * detM
        return detM
    }

    func applicationWillTerminate(_ aNotification: Notification) {
        // Insert code here to tear down your application
    }
    
    open func transpose<T>(_ input: [[T]]) -> [[T]] {
        if input.isEmpty { return [[T]]() }
        let count = input[0].count
        var out = [[T]](repeating: [T](), count: count)
        for outer in input {
            for (index, inner) in outer.enumerated() {
                out[index].append(inner)
            }
        }
        
        return out
    }


}

