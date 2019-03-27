//
//  metalTest.h
//  MatrixMult
//
//  Created by asd on 26/03/2019.
//  Copyright © 2019 Hollance. All rights reserved.
//

#ifndef metalTest_h
#define metalTest_h

#import <Foundation/Foundation.h>
#import <Accelerate/Accelerate.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

const NSUInteger N = 3000; //  complexity O(n^3), n=3000 -> x 40 metal gain

@interface MetalTest : NSObject
@property id<MTLDevice> device;
@property id<MTLCommandQueue> commandQueue;
@property NSString *name;

+(MetalTest*)init;
-(float*)testVectMultMatrixMetal;
-(float*)testVectMultMatrixCPUst;
-(float*)testVectMultMatrixCPUmt;
+(NSTimeInterval) timeIt : (void (^) (void))block;
@end

#endif /* metalTest_h */
