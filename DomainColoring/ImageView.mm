//
//  ImageView.m
//  DomainColoring
//
//  Created by asd on 01/04/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#import "ImageView.h"
#import "ImageHelper.h"
#import "ImageBuffer.h"
#import "MBEContext.h"

@implementation ImageView {
    MBEContext*ctx;
}

- (void)awakeFromNib {
    ctx = [MBEContext newContext];
}

- (void)drawRect:(NSRect)rect {
    [super drawRect:rect];
 
    int w=rect.size.width, h=rect.size.height;
    
    ImageBuffer*ibuff=[ImageBuffer initWithWidth:w Height:h];

    [ctx domainColor:ibuff.imgBuff32 w:w h:h];

    [[ibuff uiimage] drawInRect:rect];
}

@end
