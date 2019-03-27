//
//  main.m
//  VectMetal
//
//  Created by asd on 27/03/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#import <Foundation/Foundation.h>
#import "MetalTest.h"

int main(int argc, const char * argv[]) {
    @autoreleasepool {
        MetalTest *mt=[MetalTest init];
        NSTimeInterval lapm = [MetalTest timeIt: ^{ [mt testVectMultMatrixMetal]; } ];
        NSTimeInterval lapc = [MetalTest timeIt: ^{ [mt testVectMultMatrixCPUmt]; } ];

        NSLog(@"device: %@\nlap time metal: %.0f ms, cpu multithread: %.0f, ratio cpu/metal:%.2f\n",
              mt.name, lapm, lapc, lapc/lapm); // roughly x 23 in mt and 90 in st
    }
    return 0;
}
