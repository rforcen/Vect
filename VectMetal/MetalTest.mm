//
//  metalTest.m
//  MatrixMultMacOS
//
//  Created by asd on 26/03/2019.
//  Copyright © 2019 Hollance. All rights reserved.
//

#import "MetalTest.h"
#import "Vect.hpp"



const NSUInteger rowsA = N; // 5000
const NSUInteger columnsA = N; // 2000
const NSUInteger rowsB = columnsA;
const NSUInteger columnsB = N; // 3000
const NSUInteger rowsC = rowsA;
const NSUInteger columnsC = columnsB;

@implementation MetalTest {
}

+(MetalTest*)init {

    printf("n operations:%ld\n", rowsA*columnsA*columnsB);
    
    MetalTest*mt=[[super alloc]init];
    
    mt.device = MTLCreateSystemDefaultDevice();
    if (MPSSupportsMTLDevice(mt.device)) {
        mt.commandQueue = mt.device.newCommandQueue;
        mt.name=[mt.device name];
    }
    return mt;
}

-(MPSMatrixDescriptor*)vectDescritor: (Vect<float>&)v {
    return [MPSMatrixDescriptor matrixDescriptorWithDimensions: v.getRows() columns:v.getCols() rowBytes:v.colBytes() dataType:MPSDataTypeFloat32];
}
-(id<MTLBuffer>)vectBuffer: (Vect<float>&)v {
    return [_device newBufferWithBytes:v.getData() length:v.sizeBytes() options:MTLResourceOptionCPUCacheModeDefault ];
}

-(MPSMatrix*)vectMatrix: (Vect<float>&)v { // get descriptor & buffer and generate new matrix
    MPSMatrixDescriptor*desc=[MPSMatrixDescriptor matrixDescriptorWithDimensions: v.getRows() columns:v.getCols() rowBytes:v.colBytes() dataType:MPSDataTypeFloat32];
    id<MTLBuffer>buffer=[_device newBufferWithBytes:v.getData() length:v.sizeBytes() options:MTLResourceOptionCPUCacheModeDefault ];
    
    return [[MPSMatrix alloc] initWithBuffer: buffer descriptor:desc];
}

-(float*)testVectMultMatrixMetal {
    Vect<float>a=Vect<float>::rnd(N,N), b=Vect<float>::rnd(N,N), c=Vect<float>(N,N);
    
    id<MTLCommandBuffer>commandBuffer = [_commandQueue commandBuffer];
    MPSMatrixMultiplication *matrixMult = [ [MPSMatrixMultiplication alloc] initWithDevice: _device
                                                                 transposeLeft: false transposeRight:false
                                                                    resultRows: rowsC
                                                                 resultColumns: columnsC
                                                               interiorColumns: columnsA
                                                                         alpha: 1
                                                                          beta: 0];
    [matrixMult encodeToCommandBuffer:commandBuffer
                                      leftMatrix:[self vectMatrix:a]
                                     rightMatrix:[self vectMatrix:b]
                                    resultMatrix:[self vectMatrix:c]];
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    return  c.getData();
}

-(float*)testVectMultMatrixCPUst {
    Vect<float>a=Vect<float>::rnd(N,N), b=Vect<float>::rnd(N,N), c=a.stmatrixMult(b);
    return c.getData();
}
-(float*)testVectMultMatrixCPUmt {
    Vect<float>a=Vect<float>::rnd(N,N), b=Vect<float>::rnd(N,N), c=a.matrixMult(b);
    return c.getData();
}


-(void)getDevices {
    MTLDeviceNotificationHandler notificationHandler;
    
    //    AAPLViewController * __weak controller = self;
    notificationHandler = ^(id<MTLDevice> device, MTLDeviceNotificationName name)
    {
        //        [controller markHotPlugNotificationForDevice:device name:name];
    };
    
    id<NSObject> metalDeviceObserver = nil;
    NSArray<id<MTLDevice>> * availableDevices = MTLCopyAllDevicesWithObserver(&metalDeviceObserver,
                                                                              notificationHandler);
    printf("available METAL devices:%ld\n", [availableDevices count]);
    for (id<MTLDevice>dev in availableDevices)
        printf("device: %s", [dev.name cStringUsingEncoding:typeUTF8Text]);
}

+(NSTimeInterval) timeIt: (void (^) (void))block {
    NSDate *start = [NSDate date];
    block();
    return -[start timeIntervalSinceNow] * 1000.;
}
@end
