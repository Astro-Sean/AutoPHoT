#!/usr/bin/env python3
# -*- coding: utf-8 -*-


def run_astroscrappy(image_old,syntax):
    try:

        import time
        import logging
        import astroscrappy
        from threading import Thread

        logger = logging.getLogger(__name__)

        logger.info('Detecting/removing cosmic ray sources')

        #  is the program taking a while
        taking_while = False

        # output list
        clean_image = []

        print('Starting Astroscrappy ... ',end = '')

        # wrapper to move output to list
        def wrapper(func, args, res):

            res.append(func(*args))

        # setup astroscrappy but don't call
        def prep_astroscrappy():
            return astroscrappy.detect_cosmics(image_old.data,
                                        inmask = None,
                                        satlevel = 2**16,
                                        sepmed=False,
                                        gain = syntax['gain'],
                                        cleantype='medmask',
                                        fsmode='median')

        cray_remove_thread = Thread(target=wrapper,
                                args = (prep_astroscrappy,(),clean_image))

        # start time of astrcoscappry
        cray_time = time.time()

        # start thread
        cray_remove_thread.start()

        print('working ... ',end = '')

        # while astroscrapy is working keep alive with while loop
        while cray_remove_thread.isAlive():
            #  if it takes along time- print something to show it's not hung
            if time.time() - cray_time  > 15 and not taking_while:
                print('this may take some time ... ',end = '')
                taking_while = True

        print('done')

        clean_image = clean_image[0][1]

        return clean_image
    except Exception as e:
        logger.exception(e)
        return image_old